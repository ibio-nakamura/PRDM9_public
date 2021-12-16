import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
#matplotlibの出力をNotebook内に描画する
#%matplotlib inline
from math import log2
import sys
from collections import OrderedDict


def make_dfs_for_logo_maker(file_name):
    with open(file_name) as f:
        lines = f.readlines()
    
    all_filters_od = OrderedDict()
    current_filter = []
    pos=0
    for i, line in enumerate(lines):
        if line[0]=='>':
            if i!=0:
                all_filters_od[current_key] = current_filter
                current_key = line.lstrip('>').rstrip('\n')
                current_filter = []
                pos = 0
                continue
            else:
                current_key = line.lstrip('>').rstrip('\n')
                continue
        else:
            pos+=1
            temp_splited_list = line.replace('\n','').split('\t')
            splited_list = [float(item) for item in temp_splited_list]
            splited_list.insert(0,pos)
            current_filter.append(splited_list)

    od = OrderedDict()
    for key, value in all_filters_od.items():
        df = pd.DataFrame(value)
        df.columns = ['pos','A','C','G','T']
        df = df.set_index('pos')
        od[key] = df

    return od


def IC_logo_maker(od,key):
    PWM_df = od[key]
    PWM_df_list = PWM_df.values.tolist()
    IC_list = []

    for i, row in enumerate(PWM_df_list):
        entoropy = 0
        new_row = []
        for prob in row:
            if prob==0:
                entoropy += 0
            else:
                entoropy += -prob*log2(prob)
        entire_IC = 2-entoropy
        assert entire_IC>0, entire_IC

        for prob in row:
            new_row.append(prob*entire_IC)
        # insert position index at the first content in new row.
        new_row.insert(0, i+1)
        IC_list.append(new_row)

    IC_df = pd.DataFrame(IC_list)
    IC_df.columns = ['pos','A','C','G','T']
    IC_df = IC_df.set_index('pos')
    ww_logo = logomaker.Logo(IC_df, font_name='Arial Rounded MT Bold')
    ww_logo.style_spines(visible=False)
    ww_logo.style_spines(spines=['left', 'bottom'], visible=True)
    ww_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    ww_logo.ax.set_ylim([0,2])
    # style using Axes methods
    ww_logo.ax.set_ylabel("$I.C.$ (bits)", labelpad=-1, fontsize=15)
    ww_logo.ax.set_xlabel('position', fontsize=15)
    ww_logo.ax.xaxis.set_ticks_position('none')
    ww_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.savefig('./made_PWMs/logo_pngs/IC_{0}.png'.format(key))

    return 0


def IC_logo_maker_rc(od,key):
    PWM_df = od[key]
    PWM_df_list = PWM_df.values.tolist()
    IC_list = []
    # make PWM_df_lists reverse-complement
    rc_PWM_df_list = [[content for content in reversed(row)] for row in reversed(PWM_df_list)]

    for i, row in enumerate(rc_PWM_df_list):
        entoropy = 0
        new_row = []
        for prob in row:
            if prob==0:
                entoropy += 0
            else:
                entoropy += -prob*log2(prob)
        entire_IC = 2-entoropy
        assert entire_IC>0, entire_IC

        for prob in row:
            new_row.append(prob*entire_IC)
        new_row.insert(0, i+1)
        IC_list.append(new_row)

    IC_df = pd.DataFrame(IC_list)
    IC_df.columns = ['pos','A','C','G','T']
    IC_df = IC_df.set_index('pos')
    #print(IC_df)

    ww_logo = logomaker.Logo(IC_df, font_name='Arial Rounded MT Bold')
    ww_logo.style_spines(visible=False)
    ww_logo.style_spines(spines=['left', 'bottom'], visible=True)
    ww_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    ww_logo.ax.set_ylim([0,2])
    # style using Axes methods
    ww_logo.ax.set_ylabel("$I.C.$ (bits)", labelpad=-1, fontsize=15)
    ww_logo.ax.set_xlabel('position', fontsize=15)
    ww_logo.ax.xaxis.set_ticks_position('none')
    ww_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.savefig('./made_PWMs/logo_pngs/IC_{0}_rc.png'.format(key))

    return 0


def main(od, key):
    fig = plt.figure()
    objective_df = od[key]
    ww_logo = logomaker.Logo(objective_df, font_name='Arial Rounded MT Bold')
    ww_logo.style_spines(visible=False)
    ww_logo.style_spines(spines=['left', 'bottom'], visible=True)
    ww_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    # style using Axes methods
    ww_logo.ax.set_ylabel("weight", labelpad=-1, fontsize=15)
    ww_logo.ax.set_xlabel('position', fontsize=15)
    ww_logo.ax.xaxis.set_ticks_position('none')
    ww_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.savefig('./made_PWMs/logo_pngs/{0}.png'.format(key))
    return 0


def main_rc(od, key):
    fig = plt.figure()
    objective_df = od[key]
    new_columns = {'A':'T','C':'G','G':'C','T':'A'}
    new_index = {}
    for i in range(1,15):
        new_index[i] = 15-i
    rc_objective_df = objective_df[objective_df.columns[::-1]][::-1].rename(columns=new_columns,index=new_index)
    ww_logo = logomaker.Logo(rc_objective_df, font_name='Arial Rounded MT Bold')
    ww_logo.style_spines(visible=False)
    ww_logo.style_spines(spines=['left', 'bottom'], visible=True)
    ww_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    # style using Axes methods
    ww_logo.ax.set_ylabel("weight", labelpad=-1, fontsize=15)
    ww_logo.ax.set_xlabel('position', fontsize=15)
    ww_logo.ax.xaxis.set_ticks_position('none')
    ww_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.savefig('./made_PWMs/logo_pngs/{0}_rc.png'.format(key))
    return 0


if __name__=='__main__':
    dfs_od = make_dfs_for_logo_maker('./made_PWMs/multi_PWM.txt')
    keys = [key for key in dfs_od.keys()]
    for key in keys:
        main(od=dfs_od, key=key)
        IC_logo_maker(od=dfs_od,key=key)
        main_rc(od=dfs_od,key=key)
        IC_logo_maker_rc(od=dfs_od,key=key)
