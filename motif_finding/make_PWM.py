import sys

import json
import pandas as pd
from collections import Counter, OrderedDict
import glob
from tqdm import tqdm


def sort_for_best_dict(best_dict):
    sorted_best_dict = sorted(best_dict.items(),key=lambda x:x[1], reverse=True)
    return sorted_best_dict


def sort_for_worst_dict(worst_dict):
    sorted_worst_dict_list = sorted(worst_dict.items(),key=lambda x:x[1])
    return sorted_worst_dict_list


def make_PWM(sorted_list,top_percentage):
    represent_seq = sorted_list[0]
    selected_list = sorted_list[:int((len(sorted_list)*(top_percentage/100)))]
    selected_seqs = [[base for base in seq] for seq in selected_list]
    df = pd.DataFrame(selected_seqs)
    print(df)
    rows = []
    for pos in range(len(represent_seq)):
        column = df[pos]
        c = Counter(column)
        # add pseudo-count.
        A, C, G, T = c['A']+0.25, c['C']+0.25, c['G']+0.25, c['T']+0.25
        total = A+C+G+T
        row = [A/total,C/total,G/total,T/total]
        rows.append(row)
    
    return rows


def write_PWM_to_file(PWM_rows):
    with open('./made_PWMs/test_PWM.txt','w') as f:
        string = '>best_PWM_by_filter0\n'
        for row in PWM_rows:
            for i, item in enumerate(row):
                if i==3:
                    item = str(item)+'\n'
                else:
                    item = str(item)+'\t'
                string += item
        f.write(string)

    return 0


def write_multi_PWM_to_file(PWM_rows_dict):
    string = ''
    keys = [key for key in PWM_rows_dict]
    with open('./made_PWMs/multi_PWM.txt','w') as f:
        for key in keys:
            string += '>{0}\n'.format(key)
            PWM_rows = PWM_rows_dict[key]
            for row in PWM_rows:
                for i, item in enumerate(row):
                    if i==3:
                        item = str(item)+'\n'
                    else:
                        item = str(item)+'\t'
                    string += item
        string += '>end'
        f.write(string)

    return 0


def make_redundant_list_from_non_redundant_dict(sorted_dict, filter_num, which):
    ret_seq_list = []
    dir = './best_and_worst_dicts/count_dict/'
    if which=='best':
        count_dict = json.load(open(dir+'best_counts_filter{0}.json'.format(filter_num)))
    elif which=='worst':
        count_dict = json.load(open(dir+'worst_counts_filter{0}.json'.format(filter_num)))
    for item in sorted_dict:
        seq = item[0]
        count = count_dict[seq]
        for _ in range(count):
            ret_seq_list.append(seq)

    return ret_seq_list


def main():
    dir = './best_and_worst_dicts/nonredundant_dict/'
    ret_od = OrderedDict()
    for i in tqdm(range(88)):
        best_fname = 'best_seqs_filter{0}.json'.format(i)
        sorted_best_dict_list = sort_for_best_dict(json.load(open(dir+best_fname)))
        redundant_best_seq_list = make_redundant_list_from_non_redundant_dict(sorted_dict=sorted_best_dict_list, filter_num=i, which='best')
        ret_od['best_{0}'.format(i)] = make_PWM(sorted_list=redundant_best_seq_list,top_percentage=0.1)
    write_multi_PWM_to_file(PWM_rows_dict=ret_od)


if __name__=='__main__':
    main()
