import json
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
import math
import matplotlib.pyplot as plt
from scipy import stats
import sys
from collections import OrderedDict


def ave(rec_list):
    if len(rec_list) == 0:
        return np.nan
    else:
        return sum(rec_list)/len(rec_list)


def Noverlap_hot_cold_dict_maker(tsv_fname, hot_cold_judged_fname, chr_index):
    od = OrderedDict()
    # opening bed file which has information of start pos and prediction probability.
    df = pd.read_csv(tsv_fname, delimiter='\t')
    start_pos_list = list(df['start_pos'])
    probability_score_list = list(df['mc_prediction_probability'])
    assert len(start_pos_list) == len(probability_score_list)
    # opening json file which has information each fragment has a peak(1) or not(0).
    json_open = open(hot_cold_judged_fname)
    hot_cold_judged_dict = json.load(json_open)
    
    for i, pos in enumerate(start_pos_list):
        judgement = hot_cold_judged_dict[str(pos)]
        if judgement==1:
            key = ('chr'+chr_index, str(pos), 'hot')
            od[key] = probability_score_list[i]
        else:
            key = ('chr'+chr_index, str(pos), 'cold')
            od[key] = probability_score_list[i]

    return od


def Noverlap_PWM_hot_cold_dict_maker(PWM_score_fname, hot_cold_judged_fname, chr_index):
    od = OrderedDict()
    # opening json file which has PWM max-score.
    json_open = open(PWM_score_fname)
    json_dict = json.load(json_open)
    # opening json file which has information each fragment has a peak(1) or not(0).
    json_open = open(hot_cold_judged_fname)
    hot_cold_judged_dict = json.load(json_open)

    # fragment index is like "Chr1FragmentStartsFrom10234 < unknown description"
    for fragment_index, score in json_dict.items():
        splited_fragment_index = fragment_index.split()[0].split('FragmentStartsFrom')
        index_chr_index = splited_fragment_index[0].replace('Chr','')
        index_start_pos = splited_fragment_index[1]
        assert index_chr_index == chr_index, str(index_chr_index)+'≠'+chr_index
        
        judgement = hot_cold_judged_dict[str(index_start_pos)]
        if judgement==1:
            key = ('chr'+chr_index, str(index_start_pos), 'hot')
        else:
            key = ('chr'+chr_index, str(index_start_pos), 'cold')

        od[key] = score
    return od


def Noverlap_make_rec_chr_dict(path,chr_index):
    fname = 'rec_rate_100bp_chr{0}.json'.format(chr_index)
    jopen = open(path+fname)
    rec_dict = json.load(jopen)
    print('N of fragments ', len(rec_dict))

    # rec_chr_dict is modified rec_dict has format of {(pos, chr_index):rec_rate}
    rec_chr_dict = OrderedDict()
    for pos, rec_rate in rec_dict.items():
        new_key = (pos, chr_index)
        rec_chr_dict[new_key] = rec_rate

    return rec_chr_dict


#{1st percentile:[rec, rec,..],2st percentile:[rec,..]...}みたいな辞書を作ってくれる
def sort_and_make_lists(all_score_dict, all_rec_dict):
    info_list = []
    rec_list = []
    score_list = []
    # dictionary is sorted by scores.
    score_sorted = sorted(all_score_dict.items(), key=lambda x:x[1])

    for kvtupple in reversed(score_sorted):
        # key has 3 items within its tupple.; OrderedDict format is like [((chr1, 301, 'hot'),0.999)...]
        key = kvtupple[0]
        key_chr = str(key[0]).replace('chr','')
        key_pos = key[1]
        key_which_seq = str(key[2])
        score_value = kvtupple[1]
        
        key4rec_list = (key_pos, key_chr)
        try:
            rec_list.append(all_rec_dict[key4rec_list])
            info_list.append(key)
            score_list.append(score_value)
        except:
            continue
    assert len(info_list)==len(rec_list)
    assert len(info_list)==len(score_list)
    
    return [info_list, rec_list, score_list]