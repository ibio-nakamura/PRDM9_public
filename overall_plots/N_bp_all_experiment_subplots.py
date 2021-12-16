import sys
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from tqdm import tqdm
import math
import json
from collections import OrderedDict

from N_bp_Hot_and_Cold import Hot_and_Cold_main, PWM_Hot_and_Cold_main
from box_plot_modules import PWM_dict_maker, DL_dict_maker, \
                             make_percentile_dict, make_percent_dict
from N_bp_make_2by2_table_modules import Noverlap_make_rec_chr_dict


def ave(rec_list):
    if len(rec_list) == 0:
        return np.nan
    else:
        return math.fsum(rec_list)/len(rec_list)


def make_rec_chr_dict(chr_index):
    # rec_dict has format of {pos:rec.rate}.
    jopen = open('../calc_rec_rate_starts_from_0/rec_rate_100bp_starts_from_0/rec_rate_100bp_chr{0}.json'.format(chr_index))
    rec_dict = json.load(jopen)
    print("N of fragments: ", len(rec_dict))

    # rec_chr_dict is modified rec_dict has format of {(pos, chr_index):rec_rate}
    rec_chr_dict = {}
    for pos, rec_rate in rec_dict.items():
        new_key = (pos, chr_index)
        rec_chr_dict[new_key] = rec_rate

    return rec_chr_dict


def ret_percentile_dict(all_hot_average,all_cold_average,CNN_data_path,PWM_data_path,chr_indices,all_chr_rec_dict):
    all_DL_score_dict = {}
    all_PWM_score_dict = {}
    for chr_index in chr_indices:
        CNN_fname = 'CNN_prediction_Noverlap_3aug_chr_{0}.bed'.format(chr_index)
        PWM_fname = 'integrated_Noverlap_chr{0}.json'.format(chr_index)
        current_DL_score_dict = DL_dict_maker(CNN_data_path+CNN_fname, chr_index)
        current_PWM_score_dict = PWM_dict_maker(PWM_data_path+PWM_fname, chr_index)

        for key, value in current_DL_score_dict.items():
            all_DL_score_dict[key] = value
        for key, value in current_PWM_score_dict.items():
            all_PWM_score_dict[key] = value

    print('Length of all_DL_score_dict: {0}'.format(len(all_DL_score_dict)))
    print('Length of all_PWM_score_dict: {0}'.format(len(all_PWM_score_dict)))
    # each dictionary is going to be sorted and be percentile_dict.
    PWM_percentile_dict = make_percentile_dict(all_PWM_score_dict, all_chr_rec_dict)
    DL_percentile_dict = make_percentile_dict(all_DL_score_dict, all_chr_rec_dict)
    PWM_percent_dict = make_percent_dict(percentile_dict=PWM_percentile_dict, percentile_index=1)
    DL_percent_dict = make_percent_dict(percentile_dict=DL_percentile_dict, percentile_index=1)
    


    return [DL_percentile_dict, PWM_percentile_dict, DL_percent_dict, PWM_percent_dict]


def get_standard_error(data_array):
    N = len(data_array)
    mean = np.average(data_array)
    variance = sum([(item-mean)**2 for item in data_array])/(N-1)
    std = np.sqrt(variance)
    se = std/np.sqrt(N)
    return se


def main():
    chr_indices = [str(i) for i in range(1,23)]
    
    print('Now making all_chr_rec_dict..')
    all_chr_rec_dict = OrderedDict()
    for chr_index in chr_indices:
        current_chr_rec_dict = Noverlap_make_rec_chr_dict(path='../calc_rec_rate_starts_from_0/rec_rate_100bp_starts_from_0/',\
                                                          chr_index=chr_index)
        #-> returned OrderedDict format has format of {(pos, chr_index):rec_rate}
        for key, value in current_chr_rec_dict.items():
            all_chr_rec_dict[key] = value

    #https://bellcurve.jp/statistics/course/8616.htmlより
    cold_percentile_dict, cold_percent_dict, all_hot_average, all_cold_average \
                                            = Hot_and_Cold_main(path_to_CNN_dir='../CNN_dir/',\
                                                                all_chr_rec_od=all_chr_rec_dict,\
                                                                chr_indices=chr_indices)

    PWM_cold_percentile_dict, PWM_cold_percent_dict, all_hot_average_2, all_cold_average_2 \
                                            = PWM_Hot_and_Cold_main(path_to_CNN_dir='../CNN_dir/',\
                                                                    all_chr_rec_od=all_chr_rec_dict,\
                                                                    chr_indices=chr_indices)

    CNN_percentile_dict, PWM_percentile_dict, CNN_percent_dict, PWM_percent_dict \
                                            = ret_percentile_dict(all_hot_average=all_hot_average,\
                                                                  all_cold_average=all_cold_average,\
                                                                  CNN_data_path='../CNN_dir//Noverlap_predictions_bed/',\
                                                                  PWM_data_path='../PWM_dir/global_comp_each_Noverlap_calc_at_34/out_jsons/',\
                                                                  chr_indices=chr_indices,\
                                                                  all_chr_rec_dict=all_chr_rec_dict)

    assert all_hot_average==all_hot_average_2, '{0}!={1}'.format(all_hot_average,all_hot_average_2)
    assert all_cold_average==all_cold_average_2

    percent_keys = [str(1*(i_percent))+("~")+str(1*(i_percent+1))+"%" for i_percent in range(10)]
    percentile_keys = [str(10*(i_percentile))+("~")+str(10*(i_percentile+1))+"%" for i_percentile in range(10)]
    CNN_percent_values = [CNN_percent_dict[key] for key in percent_keys]
    CNN_percentile_values = [CNN_percentile_dict[key] for key in percentile_keys]
    PWM_percent_values = [PWM_percent_dict[key] for key in percent_keys]
    PWM_percentile_values = [PWM_percentile_dict[key] for key in percentile_keys]
    cold_percent_values = [cold_percent_dict[key] for key in percent_keys]
    cold_percentile_values = [cold_percentile_dict[key] for key in percentile_keys]
    PWM_cold_percent_values = [PWM_cold_percent_dict[key] for key in percent_keys]
    PWM_cold_percentile_values = [PWM_cold_percentile_dict[key] for key in percentile_keys]
    
    ave_od = OrderedDict()
    ave_od['hot'] = all_hot_average
    ave_od['cold'] = all_cold_average

    CNN_od = OrderedDict()
    CNN_od['percent'] = CNN_percent_values
    CNN_od['percentile'] = CNN_percentile_values
    CNN_od['cold_percent'] = cold_percent_values
    CNN_od['cold_percentile'] = cold_percentile_values

    PWM_od = OrderedDict()
    PWM_od['percent'] = PWM_percent_values
    PWM_od['percentile'] = PWM_percentile_values
    PWM_od['cold_percent'] = PWM_cold_percent_values
    PWM_od['cold_percentile'] = PWM_cold_percentile_values

    with open('./plot_ods/average_plot_od.json','w') as f_ave:
        json.dump(ave_od,f_ave)
    with open('./plot_ods/CNN_plot_od.json','w') as f_CNN:
        json.dump(CNN_od,f_CNN)
    with open('./plot_ods/PWM_plot_od.json','w') as f_PWM:
        json.dump(PWM_od,f_PWM)


    x = np.arange(10)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax1.plot(x, [all_hot_average]*10, linestyle='--', linewidth=1,color='red', label='average of all ChIP-seq-posive')
    ax1.plot(x, [all_cold_average]*10, linestyle='--', linewidth=1 ,color='blue',label='average of all ChIP-seq-negative')
    ax1.errorbar(x, [np.average(a_list) for a_list in CNN_percent_values],\
                 yerr=[get_standard_error(a_list) for a_list in CNN_percent_values],\
                 label='CNN', linewidth=1, color='orange') 
    ax1.errorbar(x, [np.average(a_list) for a_list in PWM_percent_values],\
                 yerr=[get_standard_error(a_list) for a_list in PWM_percent_values],\
                 label='PWM',  linewidth=1, color='green')
    ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    plt.ylabel('Average recombination rate', fontsize=10)
    plt.xlabel('Prediction scores (descending order)',fontsize=10)
    plt.ylim(0,3.55*10**(-6))
    plt.xlim(-0.3,9.3)
    plt.xticks(x,["0~1%", "", "", "", "4~5%", "", "", "", "", "9~10%"])
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=10)
    #plt.legend(bbox_to_anchor=(1,0.91), loc='upper right', borderaxespad=0.5, fontsize=7)

    ax2 = fig.add_subplot(1,2,2)
    ax2.plot(x, [all_hot_average]*10, linestyle='--', linewidth=1,color='red', label='average of all ChIP-seq-posive')
    ax2.plot(x, [all_cold_average]*10, linestyle='--', linewidth=1 ,color='blue',label='average of all ChIP-seq-negative')
    ax2.errorbar(x, [np.average(a_list) for a_list in CNN_percentile_values],\
                 yerr=[get_standard_error(a_list) for a_list in CNN_percentile_values],\
                 label='CNN', linewidth=1, color='orange')
    ax2.errorbar(x, [np.average(a_list) for a_list in PWM_percentile_values],\
                 yerr=[get_standard_error(a_list) for a_list in PWM_percentile_values],\
                 label='PWM', linewidth=1, color='green')
    ax2.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    plt.xlabel('Prediction scores (descending order)', fontsize=10)
    plt.ylim(0,3.55*10**(-6))
    plt.xlim(-0.3,9.3)
    plt.xticks(x,["0~10%", "", "", "", "40~50%", "", "", "", "", "90~100%"])
    #plt.legend(bbox_to_anchor=(1,0.91), loc='upper right', borderaxespad=0.5, fontsize=7)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=10)
    plt.savefig("./100bp_png/all_autosome_subplot_starts_from_0.png")


    x = np.arange(10)
    fig = plt.figure()
    ax3 = fig.add_subplot(1,2,1)
    ax3.plot(x, [all_hot_average]*10, linestyle='--', linewidth=1,color='red', label='average of all ChIP-seq-posive')
    ax3.plot(x, [all_cold_average]*10, linestyle='--', linewidth=1 ,color='blue',label='average of all ChIP-seq-negative')

    ax3.errorbar(x, [np.average(a_list) for a_list in cold_percent_values],\
                 yerr=[get_standard_error(a_list) for a_list in cold_percent_values],\
                 label='CNN in all ChIP-seq-negative', linewidth=1, color='purple')
    ax3.errorbar(x, [np.average(a_list) for a_list in PWM_cold_percent_values],\
                 yerr=[get_standard_error(a_list) for a_list in PWM_cold_percent_values],\
                 label='PWM in all ChIP-seq-negative', linewidth=1)
    plt.xlabel('Prediction scores (descending order)',fontsize=10)
    plt.ylabel('Average recombination rate',fontsize=10)
    plt.ylim(0,3.55*10**(-6))
    plt.xlim(-0.3,9.3)
    #plt.xticks(x,["0~10%", "10~20%", "20~30%", "30~40%", "40~50%", "50~60%", "60~70%", "70~80%", "80~90%", "90~100%"])
    plt.xticks(x,["0~1%", "", "", "", "4~5%", "", "", "", "", "9~10%"])
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=10)
    #plt.legend(bbox_to_anchor=(1,0.91), loc='upper right', borderaxespad=0.5, fontsize=7)

    ax4 = fig.add_subplot(1,2,2)
    ax4.plot(x, [all_hot_average]*10, linestyle='--', linewidth=1,color='red', label='average of all ChIP-seq-posive')
    ax4.plot(x, [all_cold_average]*10, linestyle='--', linewidth=1 ,color='blue',label='average of all ChIP-seq-negative')
    ax4.errorbar(x, [np.average(a_list) for a_list in cold_percentile_values],\
                 yerr=[get_standard_error(a_list) for a_list in cold_percentile_values],\
                 label='CNN in all ChIP-seq-negative', linewidth=1, color='purple')
    ax4.errorbar(x, [np.average(a_list) for a_list in PWM_cold_percentile_values],\
                 yerr=[get_standard_error(a_list) for a_list in PWM_cold_percentile_values],\
                 label='PWM in all ChIP-seq-negative', linewidth=1)
    plt.xlabel('Prediction scores (descending order)', fontsize=10)
    plt.ylim(0,3.55*10**(-6))
    plt.xlim(-0.3,9.3)
    plt.xticks(x,["0~10%", "", "", "", "40~50%", "", "", "", "", "90~100%"])
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=10)
    #plt.legend(bbox_to_anchor=(1,0.91), loc='upper right', borderaxespad=0.5, fontsize=7)
    plt.savefig("./100bp_png/ChIP-seq-negative_autosome_subplot_starts_from_0.png")

    print(len(CNN_percentile_values[0]))
    print(ave(CNN_percentile_values[0]))
    return 0


if __name__=='__main__':
    main()