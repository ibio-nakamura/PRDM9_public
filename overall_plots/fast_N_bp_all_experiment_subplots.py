import json
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from tqdm import tqdm
import math
from collections import OrderedDict


def get_standard_error(data_array):
    N = len(data_array)
    mean = np.average(data_array)
    #print('mean: ', mean)
    variance = sum([(item-mean)**2 for item in data_array])/(N-1)
    std = np.sqrt(variance)
    #print('std: ', std)
    se = std/np.sqrt(N)
    return se


CNN_od = json.load(open('./plot_ods/CNN_plot_od.json'))
PWM_od = json.load(open('./plot_ods/PWM_plot_od.json'))
ave_od = json.load(open('./plot_ods/average_plot_od.json'))
all_hot_average = ave_od['hot']
all_cold_average = ave_od['cold']
CNN_percent_values = [i for i in CNN_od['percent']]
CNN_percentile_values = [j for j in CNN_od['percentile']]
CNN_cold_percent_values = [k for k in CNN_od['cold_percent']]
CNN_cold_percentile_values = [l for l in CNN_od['cold_percentile']]
PWM_percent_values = [i for i in PWM_od['percent']]
PWM_percentile_values = [j for j in PWM_od['percentile']]
PWM_cold_percent_values = [k for k in PWM_od['cold_percent']]
PWM_cold_percentile_values = [l for l in PWM_od['cold_percentile']]

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
plt.ylabel('Average recombination rate (cM/bp)', fontsize=10)
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
plt.savefig("./plot_ods/fast_all_autosome_subplot_starts_from_0.png")
x = np.arange(10)
fig = plt.figure()
ax3 = fig.add_subplot(1,2,1)
ax3.plot(x, [all_hot_average]*10, linestyle='--', linewidth=1,color='red', label='average of all ChIP-seq-posive')
ax3.plot(x, [all_cold_average]*10, linestyle='--', linewidth=1 ,color='blue',label='average of all ChIP-seq-negative')
ax3.errorbar(x, [np.average(a_list) for a_list in CNN_cold_percent_values],\
             yerr=[get_standard_error(a_list) for a_list in CNN_cold_percent_values],\
             label='CNN in all ChIP-seq-negative', linewidth=1, color='purple')
ax3.errorbar(x, [np.average(a_list) for a_list in PWM_cold_percent_values],\
             yerr=[get_standard_error(a_list) for a_list in PWM_cold_percent_values],\
             label='PWM in all ChIP-seq-negative', linewidth=1)
plt.xlabel('Prediction scores (descending order)',fontsize=10)
plt.ylabel('Average recombination rate (cM/bp)',fontsize=10)
plt.ylim(0,3.55*10**(-6))
plt.xlim(-0.3,9.3)
plt.xticks(x,["0~1%", "", "", "", "4~5%", "", "", "", "", "9~10%"])
plt.xticks(fontsize=7)
plt.yticks(fontsize=10)
#plt.legend(bbox_to_anchor=(1,0.91), loc='upper right', borderaxespad=0.5, fontsize=7)
ax4 = fig.add_subplot(1,2,2)
ax4.plot(x, [all_hot_average]*10, linestyle='--', linewidth=1,color='red', label='average of all ChIP-seq-posive')
ax4.plot(x, [all_cold_average]*10, linestyle='--', linewidth=1 ,color='blue',label='average of all ChIP-seq-negative')
ax4.errorbar(x, [np.average(a_list) for a_list in CNN_cold_percentile_values],\
             yerr=[get_standard_error(a_list) for a_list in CNN_cold_percentile_values],\
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
plt.savefig("./plot_ods/fast_ChIP-seq-negative_autosome_subplot_starts_from_0.png")