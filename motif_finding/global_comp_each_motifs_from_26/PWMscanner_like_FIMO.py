import numpy as np
import copy
from PWMmaker import PWMmaker
from data_utils import reverse_complement
from math import log2


def fragment_base_comp(fragment):
    fragment_len = len(fragment)
    fragment_base_count = [ fragment.count(base) for base in ['A','C','G','T'] ]  
    fragment_base_composition = [ (count/fragment_len) for count in fragment_base_count ]

    return fragment_base_composition


def scanner(seq, PWM):
    dict = {"A":0, "C":1, "G":2, "T":3}

    if not('N' in seq):
        likelihood_list = []
        PWM_len = len(PWM)
        seq_len = len(seq)
        #the number which shows how many times the matrix scans(note: not how many times it moves).
        scan = seq_len - PWM_len + 1
        
        for i in range(scan):
            prob_list = []
            rc_prob_list = []
            scanned_region = seq[i:i+PWM_len]

            rv_scanned_region = reverse_complement(scanned_region)
            for scanned_content in [scanned_region, rv_scanned_region]:
                for j, base in enumerate(scanned_content):
                    composition_list = fragment_base_comp(scanned_content)
                    lprob = PWM[j][dict[base]]/composition_list[dict[base]]
                    if scanned_content == scanned_region:
                        prob_list.append(lprob)
                    else:
                        rc_prob_list.append(lprob)

            log_score_list = np.log2(prob_list)
            rc_log_score_list = np.log2(rc_prob_list)
            log_likelihood = np.sum(log_score_list)
            rc_log_likelihood = np.sum(rc_log_score_list)

            bigger_log_likelihood = max(log_likelihood, rc_log_likelihood)
            likelihood_list.append(bigger_log_likelihood)
        return max(likelihood_list)
    else:
        return None


def scanner_fragment_comp(seq, PWM):
    dict = {"A":0, "C":1, "G":2, "T":3}
    rc_seq = reverse_complement(seq)
    seq_comp = fragment_base_comp(seq)
    rc_seq_comp = fragment_base_comp(rc_seq)

    if not('N' in seq):
        likelihood_list = []
        PWM_len = len(PWM)
        seq_len = len(seq)
        #the number which shows how many times the matrix scans(note: not how many times it moves).
        scan = seq_len - PWM_len + 1
        
        for i in range(scan):
            prob_list = []
            rc_prob_list = []
            scanned_region = seq[i:i+PWM_len]

            rv_scanned_region = reverse_complement(scanned_region)
            for scanned_content in [scanned_region, rv_scanned_region]:
                for j, base in enumerate(scanned_content):
                    if scanned_content == scanned_region:
                        lprob = PWM[j][dict[base]]/seq_comp[dict[base]]
                        prob_list.append(lprob)
                    else:
                        lprob = PWM[j][dict[base]]/rc_seq_comp[dict[base]]
                        rc_prob_list.append(lprob)

            log_score_list = np.log2(prob_list)
            rc_log_score_list = np.log2(rc_prob_list)
            log_likelihood = np.sum(log_score_list)
            rc_log_likelihood = np.sum(rc_log_score_list)

            bigger_log_likelihood = max(log_likelihood, rc_log_likelihood)
            likelihood_list.append(bigger_log_likelihood)
        return max(likelihood_list)
    else:
        return None


def scanner_can_get_base_composition(seq, PWM, base_composition_list, rc_base_composition_list):
    dict = {"A":0, "C":1, "G":2, "T":3}

    if not('N' in seq):
        likelihood_list = []
        PWM_len = len(PWM)
        seq_len = len(seq)
        #the number which shows how many times the matrix scans(note: not how many times it moves).
        scan = seq_len - PWM_len + 1
        
        for i in range(scan):
            prob_list = []
            rc_prob_list = []
            scanned_region = seq[i:i+PWM_len]
            rv_scanned_region = reverse_complement(scanned_region)
            for scanned_content in [scanned_region, rv_scanned_region]:
                for j, base in enumerate(scanned_content):
                    if scanned_content == scanned_region:
                        lprob = PWM[j][dict[base]]/base_composition_list[dict[base]]
                        if lprob==0:
                            break
                        else:
                            prob_list.append(lprob)
                    else:
                        lprob = PWM[j][dict[base]]/rc_base_composition_list[dict[base]]
                        if lprob==0:
                            break
                        else:
                            rc_prob_list.append(lprob)

            log_score_list = np.log2(prob_list)
            rc_log_score_list = np.log2(rc_prob_list)
            log_likelihood = np.sum(log_score_list)
            rc_log_likelihood = np.sum(rc_log_score_list)

            bigger_log_likelihood = max(log_likelihood, rc_log_likelihood)
            likelihood_list.append(bigger_log_likelihood)
        return max(likelihood_list)
    else:
        return None
