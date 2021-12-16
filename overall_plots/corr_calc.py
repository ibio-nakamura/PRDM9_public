import sys

import json
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
import math
from scipy import stats
from collections import OrderedDict

from box_plot_modules import PWM_dict_maker, DL_dict_maker


def make_rec_chr_dict(chr_index):
    # rec_dict has format of {pos:rec.rate}.
    fname = 'rec_rate_100bp_chr{0}.json'.format(chr_index)
    jopen = open('../calc_rec_rate_starts_from_0/rec_rate_100bp_starts_from_0/'+fname)
    rec_dict = json.load(jopen)
    print("N of fragments: ", len(rec_dict))

    # rec_chr_dict is modified rec_dict has format of {(pos, chr_index):rec_rate}
    rec_chr_dict = OrderedDict()
    for pos, rec_rate in rec_dict.items():
        new_key = (pos, chr_index)
        rec_chr_dict[new_key] = rec_rate

    return rec_chr_dict


def score_and_rec_sort(score_dict, rec_dict):
    score_sorted = sorted(score_dict.items(), key=lambda x:x[1])
    rec_list = []
    score_list = []
    for kvtupple in reversed(score_sorted):
        key = kvtupple[0]
        score_value = kvtupple[1]
        
        try:
            rec_list.append(rec_dict[key])
            score_list.append(score_value)
        except KeyError:
            continue
    
    return [score_list, rec_list]


def main4all():
    all_rec_dict = {}
    all_CNN_score_dict = {}
    all_PWM_score_dict = {}
    chr_indices = [str(i) for i in range(1,23)]
    for chr_index in chr_indices:
        CNN_path = '../CNN_dir/Noverlap_predictions_bed/'
        CNN_fname = 'CNN_prediction_Noverlap_3aug_chr_{0}.bed'.format(chr_index)
        PWM_path = '../PWM_dir/global_comp_each_Noverlap_calc_at_34/out_jsons/'
        PWM_fname = 'integrated_Noverlap_chr{0}.json'.format(chr_index)

        current_CNN_score_dict = DL_dict_maker(CNN_path+CNN_fname, chr_index)
        current_PWM_score_dict = PWM_dict_maker(PWM_path+PWM_fname, chr_index)
        current_rec_dict = make_rec_chr_dict(chr_index)

        for key, value in current_rec_dict.items():
            all_rec_dict[key] = value
        for key, value in current_CNN_score_dict.items():
            all_CNN_score_dict[key] = value
        for key, value in current_PWM_score_dict.items():
            all_PWM_score_dict[key] = value

    CNN_score_list, CNN_rec_list = score_and_rec_sort(score_dict=all_CNN_score_dict, rec_dict=all_rec_dict)
    PWM_score_list, PWM_rec_list = score_and_rec_sort(score_dict=all_PWM_score_dict, rec_dict=all_rec_dict)

    assert len(CNN_score_list)==len(PWM_score_list)
    for i in range(len(CNN_score_list)-1):
        assert 0<=CNN_score_list[i]-CNN_score_list[i+1] 
        assert 0<=PWM_score_list[i]-PWM_score_list[i+1], 'sort is troubled in {0} and {1} at {2}'.format(PWM_score_list[i],CNN_score_list[i+1],i)
    
    print('CNN')
    print('spearmanr in CNN')
    print(stats.spearmanr(CNN_score_list,CNN_rec_list))
    with open('./spearmanr/CNN_all_autosomal_spearmanr.txt','w') as f:
        f.write(str(stats.spearmanr(CNN_score_list,CNN_rec_list)))

    print('PWM')
    print('spearmanr in PWM')
    print(stats.spearmanr(PWM_score_list,PWM_rec_list))
    with open('./spearmanr/PWM_all_autosomal_spearmanr.txt','w') as f:
        f.write(str(stats.spearmanr(PWM_score_list,PWM_rec_list)))

    print('data: ',len(PWM_score_list))


if __name__=='__main__':
    main4all()
