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


def PWM_dict_maker(fname, chr_index):
    ret_od = OrderedDict()
    jopen = open(fname)
    original_dict = json.load(jopen)
    for start_key, value in original_dict.items():
        #print('start_key: {0}'.format(start_key))
        chr_part, pos_part = start_key.split('FragmentStartsFrom')
        chr_part = chr_part.replace('Chr','')
        pos_part = pos_part.replace(' <unknown description>','')
        ret_od[(pos_part,chr_part)] = value

    return ret_od


#DLはbed形式でスコアが保存されている。
def DL_dict_maker(DL_fname, chr_index):
    ret_dict = OrderedDict()
    df = pd.read_csv(DL_fname, delimiter='\t')
    start_pos_list = list(df['start_pos'])
    score_list = list(df['mc_prediction_probability'])
    for start, score in zip(start_pos_list, score_list):
        ret_dict[(str(start),chr_index)] = score
    
    return ret_dict


def make_percentile_dict(score_dict, all_rec_dict):
    ret_dict = OrderedDict()
    score_sorted = sorted(score_dict.items(), key=lambda x:x[1])
    rec_list = []
    score_list = []
    for kvtupple in reversed(score_sorted):
        key = kvtupple[0]
        score_value = kvtupple[1]
        try:
            rec_list.append(all_rec_dict[key])
            score_list.append(score_value)
        except KeyError:
            continue

    print("average whole recombination rate: ", (np.mean(rec_list)))
    devided_rec_list = list(np.array_split(rec_list,10))
    devided_score_list = list(np.array_split(score_list,10))
    devided_rec_list = [list(one_percentile) for one_percentile in devided_rec_list]
    devided_score_list = [list(one_percentile) for one_percentile in devided_score_list] 

    for i_percentile, rec_list_in_a_percentile in enumerate(devided_rec_list):
        ret_dict[str(10*(i_percentile))+("~")+str(10*(i_percentile+1))+"%"] = rec_list_in_a_percentile
        #print('length of {0}-th percentile is {1}'.format(i_percentile,len(rec_list_in_a_percentile)))

    od = OrderedDict()
    od['percentile_score_list'] = devided_score_list
    od['percentile_rec_list'] = devided_rec_list
    with open('./record/CNN_percentile_dict.json', 'w') as f:
        json.dump(od,f)

    return ret_dict


#上のmake_percentile_dictを受けて、
#{"0~1%":[rec, rec,..],"1~2%":[rec,..]...}みたいな辞書を作ってくれる
def make_percent_dict(percentile_dict,percentile_index=1):
    ret_dict = OrderedDict()
    rec_key = str(10*(percentile_index-1))+"~"+str(10*(percentile_index))+"%"
    print('Now making percent_dict with key of {0} in percentile_dict'.format(rec_key))
    rec_list_in_a_percentile = percentile_dict[rec_key]
    devided_rec_list = list(np.array_split(rec_list_in_a_percentile,10))
    devided_rec_list = [list(one_percentile) for one_percentile in devided_rec_list]
    for i_percent, rec_list_in_a_percent in enumerate(devided_rec_list):
        ret_dict[str(1*(i_percent))+("~")+str(1*(i_percent+1))+"%"]=rec_list_in_a_percent
        #print('length of {0}-th percent is {1}'.format(i_percent,len(rec_list_in_a_percent)))
    
    od = OrderedDict()
    od['percent_rec_list'] = devided_rec_list
    with open('./record/CNN_percent_dict.json', 'w') as f:
        json.dump(od,f)

    return ret_dict
