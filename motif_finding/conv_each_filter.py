import sys
import types
import json

import numpy as np
import keras
from keras.models import Sequential
from keras.initializers import Constant
import keras.backend as K
from keras.engine.topology import Layer
import tensorflow as tf
import h5py
import pandas as pd
import random
from collections import OrderedDict
from tqdm import tqdm

from data_utils import load_model_yaml, one_hot, load_fasta_gz


#reference: https://atmarkit.itmedia.co.jp/ait/articles/2005/13/news009.html
def lrelu(x, alpha=0.01):
  return np.where(x >= 0.0, x, alpha * x)


def load_input_data(input_len):
    all_sequences = []
    chr_indices = [str(i) for i in range(1,23)]
    chr_indices.append('X')
    for chr_index in chr_indices:
        PATH = '../make_Noverlap_fragments/100bp_Noverlap_fastas/'
        current_sequences = load_fasta_gz(f_name=PATH+'100bp_NoverlapChr{0}.fasta.gz'.format(chr_index), \
                                          input_len=input_len)
        current_sequences = [seq for seq in current_sequences if not('N' in seq)]
        all_sequences += current_sequences
    
    sample_size = len(all_sequences)//100
    sampled_sequences = random.sample(all_sequences,sample_size)
    print('N of sequences: ', len(sampled_sequences))
    onehot_samples = np.array([one_hot(seq) for seq in sampled_sequences])

    return onehot_samples


def handmade_conv(input,filter_num,filter_length,position):
    with h5py.File('../CNN_dir/SavedModel_100bp_3aug.h5', 'r') as h5:
        kernel = h5["conv1d_1/conv1d_1/kernel:0"]
        biases = h5["conv1d_1/conv1d_1/bias:0"]
        bias = np.array(biases)[filter_num]
        filter = np.array(kernel)[:,:,filter_num]
    
    conved_range_seq = input[position:filter_length+position]
    sum = np.sum(conved_range_seq*filter) + bias
    score_before_activation = sum
    #score = lrelu(sum)
    #print("output of filter of position{0}".format(position))
    
    return [score_before_activation,conved_range_seq]


def handmade_conv_automatic_scan(input,filter,bias):
    filter_length = len(filter)
    input_len = len(input)
    for position in range(input_len-filter_length+1):
        conved_range_seq = input[position:filter_length+position]
        sum = np.sum(conved_range_seq*filter) + bias
        score_before_activation = sum
        if position==0:
            best_score = score_before_activation
            best_sequence = conved_range_seq
            worst_score = score_before_activation
            worst_sequence = conved_range_seq
        if best_score <= score_before_activation:
            best_score = score_before_activation
            best_sequence = conved_range_seq
        if score_before_activation <= worst_score:
            worst_score = score_before_activation
            worst_sequence = conved_range_seq

        #score = lrelu(sum)
        #print("output of filter of position{0}".format(position))
    
    return [worst_score,worst_sequence,best_score,best_sequence]


def arc_one_hot(array):
    seq = ''
    """
    one_hot_conv = {"A": [1,0,0,0],
                    "T": [0,0,0,1],
                    "C": [0,1,0,0],
                    "G": [0,0,1,0]}
    """
    arc_one_hot_conv = {0:'A',
                        3:'T',
                        1:'C',
                        2:'G'}
    string_list = [arc_one_hot_conv[np.where(one_hot_vec==1)[0][0]] for one_hot_vec in array]
    for string in string_list:
        seq += string

    return seq


if __name__=='__main__':
    input_len = 100
    seqs = load_input_data(input_len)
    best_params_dict = json.load(open('../CNN_dir/best_params_3aug.json'))
    all_filter_num = best_params_dict['filter_num']
    filter_length = best_params_dict['first_filter_length']
    for filter_num in range(all_filter_num):
        best_score_and_seq_od = OrderedDict()
        worst_score_and_seq_od = OrderedDict()
        best_count_dict = OrderedDict()
        worst_count_dict = OrderedDict()
        with h5py.File('../CNN_dir/SavedModel_100bp_3aug.h5', 'r') as h5:
            kernel = h5["conv1d_1/conv1d_1/kernel:0"]
            biases = h5["conv1d_1/conv1d_1/bias:0"]
            bias = np.array(biases)[filter_num]
            filter = np.array(kernel)[:,:,filter_num]
            print(filter)
        for input in tqdm(seqs):
            worst_score,worst_sequence,best_score,best_sequence = handmade_conv_automatic_scan(input,filter,bias)
            if arc_one_hot(best_sequence) in best_score_and_seq_od:
                assert best_score_and_seq_od[arc_one_hot(best_sequence)] == float(best_score)
                best_count_dict[arc_one_hot(best_sequence)] += 1
            else:
                best_score_and_seq_od[arc_one_hot(best_sequence)] = float(best_score)
                best_count_dict[arc_one_hot(best_sequence)] = 1

            if arc_one_hot(worst_sequence) in worst_score_and_seq_od:
                assert worst_score_and_seq_od[arc_one_hot(worst_sequence)] == float(worst_score)
                worst_count_dict[arc_one_hot(worst_sequence)] += 1
            else:
                worst_score_and_seq_od[arc_one_hot(worst_sequence)] = float(worst_score)
                worst_count_dict[arc_one_hot(worst_sequence)] = 1
        
        # nonredundant_dict has scores of sequence.
        with open('./best_and_worst_dicts/nonredundant_dict/best_seqs_filter{0}.json'.format(filter_num),'w') as f:
            json.dump(best_score_and_seq_od, f, indent=4)
        with open('./best_and_worst_dicts/nonredundant_dict/worst_seqs_filter{0}.json'.format(filter_num),'w') as f:
            json.dump(worst_score_and_seq_od, f, indent=4)
        # count_dict has count for each sequence.
        with open('./best_and_worst_dicts/count_dict/best_counts_filter{0}.json'.format(filter_num),'w') as f:
            json.dump(best_count_dict, f, indent=4)
        with open('./best_and_worst_dicts/count_dict/worst_counts_filter{0}.json'.format(filter_num),'w') as f:
            json.dump(worst_count_dict, f, indent=4)
