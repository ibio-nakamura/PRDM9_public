import json
from collections import Counter
from tqdm import tqdm
import numpy as np
from sklearn.metrics import roc_curve
from sklearn import metrics


def return_AUC(y_true,preds):
    fpr, tpr, thresholds = roc_curve(y_true, preds)
    auc = metrics.auc(fpr, tpr)

    return auc


def return_PWM_scores_from_json(chr_index,path,motif_num):
    json_fname_for_pos = "positive_chr_{0}_global_comp_made_motif{1}.json".format(chr_index,motif_num)
    json_fname_for_neg = "negative_chr_{0}_global_comp_made_motif{1}.json".format(chr_index,motif_num)
    dict_for_pos = json.load(open(path+json_fname_for_pos))
    dict_for_neg = json.load(open(path+json_fname_for_neg))
    scores_for_pos = [score for score in dict_for_pos.values()]
    scores_for_neg = [score for score in dict_for_neg.values()]

    return [scores_for_pos,scores_for_neg]


chr_indices = [str(i) for i in range(1,23)]
chr_indices.append('X')
for motif_num in range(44):
    pos_scores_for_a_motif = []
    neg_scores_for_a_motif = []
    for chr_index in chr_indices:
        current_pos_scores_for_a_motif, current_neg_scores_for_a_motif \
                = return_PWM_scores_from_json(chr_index=chr_index,path='./loglikelihood_ratio_scores/',motif_num=motif_num)
        pos_scores_for_a_motif += current_pos_scores_for_a_motif
        neg_scores_for_a_motif += current_neg_scores_for_a_motif
    
    y_true = [1]*len(pos_scores_for_a_motif)+[0]*len(neg_scores_for_a_motif)
    auc = return_AUC(y_true=y_true,preds=pos_scores_for_a_motif+neg_scores_for_a_motif)
    print('AUC of motif{0}: {1}'.format(motif_num,auc))
