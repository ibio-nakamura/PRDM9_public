import sys
import json
import math
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from scipy import stats
from collections import OrderedDict
from ipywidgets import interact, FloatSlider

from N_bp_make_2by2_table_modules import Noverlap_hot_cold_dict_maker, Noverlap_make_rec_chr_dict, \
                                         sort_and_make_lists, Noverlap_PWM_hot_cold_dict_maker
from box_plot_modules import PWM_dict_maker, DL_dict_maker, \
                             make_percentile_dict, make_percent_dict


def ave(rec_list):
    if len(rec_list) == 0:
        return np.nan
    else:
        return math.fsum(rec_list)/len(rec_list)


# すでにスコアでソートされたlistが引数になっている。
def make_sorted_hot_cold_each(CNN_info_list, CNN_rec_list, CNN_score_list):
    CNN_only_cold_rec_list = []
    CNN_only_cold_score_list = []
    CNN_only_hot_rec_list = []
    CNN_only_hot_score_list = []

    # Firstly, eliminate hot fragments-related contents in rec_list and score_list using info_list.
    for current_index, CNN_info  in enumerate(CNN_info_list):
        CNN_which_seq = CNN_info[2]
        if CNN_which_seq == 'hot':
            CNN_only_hot_rec_list.append(CNN_rec_list[current_index])
            CNN_only_hot_score_list.append(CNN_score_list[current_index])
        elif CNN_which_seq == 'cold':
            CNN_only_cold_rec_list.append(CNN_rec_list[current_index])
            CNN_only_cold_score_list.append(CNN_score_list[current_index])
    ChIP_pos_average = ave(CNN_only_hot_rec_list)
    print('ChIP fragment average from make_sorted_hot_cold_each: ', ChIP_pos_average)
    ChIP_neg_average = ave(CNN_only_cold_rec_list)
    print('ChIP negative fragment average from make_sorted_hot_cold_each: ', ChIP_neg_average)

    table = [CNN_only_hot_rec_list, CNN_only_hot_score_list, CNN_only_cold_rec_list, CNN_only_cold_score_list]
    return table


def make_percentile_dict_for_hot_cold_each(score_list, rec_list):
    assert len(score_list)==len(rec_list), '{0}!={1}'.format(len(score_list),len(rec_list))
    for i in range(len(score_list)-1):
        assert 0 <= (score_list[i] - score_list[i+1]), 'sort is troubled at {0} and {1}'.format(score_list[i],score_list[i+1])
    
    devided_rec_list = list(np.array_split(rec_list,10))
    devided_rec_list = [list(one_percentile) for one_percentile in devided_rec_list]
    ret_dict = OrderedDict()
    for i_percentile, rec_list_in_a_percentile in enumerate(devided_rec_list):
        ret_dict[str(10*(i_percentile))+("~")+str(10*(i_percentile+1))+"%"] = rec_list_in_a_percentile
    
    return ret_dict


def Hot_and_Cold_main(path_to_CNN_dir,all_chr_rec_od,chr_indices):
    CNN_all_score_od = OrderedDict()

    for chr_index in chr_indices:
        # preparing judgement file which has information a peak is on a fragment(1) or not(0).
        hot_cold_judged_fname = 'peak_judgement_out/judged_chr_{0}.json'.format(chr_index)
        # preparing DeepLearning score ordered dict...
        tsv_fname = 'Noverlap_predictions_bed/CNN_prediction_Noverlap_3aug_chr_{0}.bed'.format(chr_index)
        current_od = Noverlap_hot_cold_dict_maker(path_to_CNN_dir+tsv_fname,path_to_CNN_dir+hot_cold_judged_fname,chr_index)
        #-> returned OrderedDict format is like [((chr1, 301, 'hot'),0.999)...]
        for key, value in current_od.items():
            CNN_all_score_od[key] = value

    CNN_info_list, CNN_rec_list, CNN_score_list = sort_and_make_lists(CNN_all_score_od, all_chr_rec_od)
    od = OrderedDict()
    od['CNN_score_list']=CNN_score_list
    od['CNN_rec_list']=CNN_rec_list
    with open('./record/sorted_CNN_score.json','w') as f:
        json.dump(od,f)

    # hotとcoldをinfo_listから判定して振り分ける
    CNN_only_hot_rec_list, CNN_only_hot_score_list, CNN_only_cold_rec_list, CNN_only_cold_score_list\
                                                            = make_sorted_hot_cold_each(CNN_info_list, CNN_rec_list, CNN_score_list)
    
    assert len(CNN_only_hot_rec_list)+len(CNN_only_cold_rec_list) == len(CNN_rec_list)
    for i in range(len(CNN_only_cold_score_list)-1):
        assert 0<=CNN_only_cold_score_list[i]-CNN_only_cold_score_list[i+1],'sort is troubled at {0} and {1}'.format(CNN_only_cold_score_list[i],CNN_only_cold_score_list[i+1])
    
    ChIP_pos_percentile_dict = make_percentile_dict_for_hot_cold_each(score_list=CNN_only_hot_score_list,rec_list=CNN_only_hot_rec_list)
    ChIP_neg_percentile_dict = make_percentile_dict_for_hot_cold_each(score_list=CNN_only_cold_score_list,rec_list=CNN_only_cold_rec_list)
    all_hot_average = ave(CNN_only_hot_rec_list)
    all_cold_average = ave(CNN_only_cold_rec_list)
    print('all_hot_average: ', all_hot_average)
    print('all_cold_average: ', all_cold_average)
    cold_percent_dict = make_percent_dict(percentile_dict=ChIP_neg_percentile_dict)
    #with open('CNN_spearmanr_in_ChIP-negative.json','w') as f:
    #    f.write(str(stats.spearmanr(CNN_only_cold_rec_list, CNN_only_cold_score_list)))
    
    od = OrderedDict()
    od['cold_rec_list'] = CNN_only_cold_rec_list
    od['cold_score_list'] = CNN_only_cold_score_list
    with open('./record/CNN_spearman_data.json','w') as f:
        json.dump(od,f)

    od = OrderedDict()
    od['ChIP_pos_percentile_dict']=ChIP_pos_percentile_dict
    od['ChIP_neg_percentile_dict']=ChIP_neg_percentile_dict
    od['cold_percent_dict']=cold_percent_dict
    with open('./record/CNN_ChIP_pos_neg.json','w') as f:
        json.dump(od,f)

    return [ChIP_neg_percentile_dict, cold_percent_dict, all_hot_average, all_cold_average]


def PWM_Hot_and_Cold_main(path_to_CNN_dir,all_chr_rec_od,chr_indices):
    PWM_all_score_od = OrderedDict()

    for chr_index in chr_indices:
        # preparing judgement file which has information a peak is on a fragment(1) or not(0).
        hot_cold_judged_fname = 'peak_judgement_out/judged_chr_{0}.json'.format(chr_index)
        # preparing DeepLearning score ordered dict...
        PWM_score_fname = '../PWM_dir/global_comp_each_Noverlap_calc_at_34/out_jsons/integrated_Noverlap_chr{0}.json'.format(chr_index)
        current_od = Noverlap_PWM_hot_cold_dict_maker(PWM_score_fname, path_to_CNN_dir+hot_cold_judged_fname, chr_index)
        #-> returned OrderedDict format is like [((chr1, 301, 'hot'),0.999)...]
        for key, value in current_od.items():
            PWM_all_score_od[key] = value
    
    PWM_info_list, PWM_rec_list, PWM_score_list = sort_and_make_lists(PWM_all_score_od, all_chr_rec_od)
    od = OrderedDict()
    od['PWM_score_list']=PWM_score_list
    od['PWM_rec_list']=PWM_rec_list
    with open('./record/sorted_PWM_score.json','w') as f:
        json.dump(od,f)
    
    # Determine hot and cold from info_list and sort them.
    PWM_only_hot_rec_list, PWM_only_hot_score_list, PWM_only_cold_rec_list, PWM_only_cold_score_list\
                                                = make_sorted_hot_cold_each(PWM_info_list, PWM_rec_list, PWM_score_list)
    assert len(PWM_only_hot_rec_list)+len(PWM_only_cold_rec_list) == len(PWM_rec_list)
    for i in range(len(PWM_only_cold_score_list)-1):
        assert 0<=PWM_only_cold_score_list[i]-PWM_only_cold_score_list[i+1],'sort is troubled at {0} and {1}'.format(PWM_only_cold_score_list[i],PWM_only_cold_score_list[i+1]) 

    ChIP_pos_percentile_dict = make_percentile_dict_for_hot_cold_each(score_list=PWM_only_hot_score_list,rec_list=PWM_only_hot_rec_list)
    ChIP_neg_percentile_dict = make_percentile_dict_for_hot_cold_each(score_list=PWM_only_cold_score_list,rec_list=PWM_only_cold_rec_list)
    all_hot_average = ave(PWM_only_hot_rec_list)
    all_cold_average = ave(PWM_only_cold_rec_list)
    print('all_hot_average: ', all_hot_average)
    print('all_cold_average: ', all_cold_average)
    cold_percent_dict = make_percent_dict(percentile_dict=ChIP_neg_percentile_dict)
    #with open('PWM_spearmanr_in_ChIP-negative.json','w') as f:
    #    f.write(str(stats.spearmanr(PWM_only_cold_rec_list, PWM_only_cold_score_list)))
    od = OrderedDict()
    od['cold_rec_list'] = PWM_only_cold_rec_list
    od['cold_score_list'] = PWM_only_cold_score_list
    with open('./record/PWM_spearman_data.json','w') as f:
        json.dump(od,f)

    od = OrderedDict()
    od['ChIP_pos_percentile_dict']=ChIP_pos_percentile_dict
    od['ChIP_neg_percentile_dict']=ChIP_neg_percentile_dict
    od['cold_percent_dict']=cold_percent_dict
    with open('./record/PWM_ChIP_pos_neg.json','w') as f:
        json.dump(od,f)

    return [ChIP_neg_percentile_dict, cold_percent_dict, all_hot_average, all_cold_average]
