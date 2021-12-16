from PWMmaker import PWMmaker
from PWMscanner_like_FIMO import scanner_can_get_base_composition, fragment_base_comp
import sys
import json
from data_utils import load_fasta_gz_and_return_OrderedDict, load_fasta_gz, reverse_complement, load_hg19_gz
from tqdm import tqdm
import pandas as pd
from collections import OrderedDict


def log_likelihood_each_matrix_global(PWM_fname, chr_index, fragments_od, input_len):
    log_likelihood_od = OrderedDict()
    global_seq = load_hg19_gz(f_name='../../hg19/chr'+chr_index+'.fa.gz')[0].upper().replace('N','')
    base_composition_list = fragment_base_comp(global_seq)
    rc_base_composition_list = fragment_base_comp(reverse_complement(global_seq))

    for key in tqdm(fragments_od):
        fragment = fragments_od[key]
        PWM = PWMmaker('../made_PWMs/dev_PWM/'+PWM_fname)
        #the max log-likelihood will be returned in a fragments(301bp).
        log_likelihood = scanner_can_get_base_composition(fragment, PWM, base_composition_list, rc_base_composition_list)
        if log_likelihood != None:
            log_likelihood_od[key] = log_likelihood

    return log_likelihood_od


# In main function, you can make json file with expectation of the number of motifs and position taking account of fragment has 'N'.
def main(chr_index,PWM_num,input_len):
    PWM_fname = 'best_'+PWM_num+'.txt'
    for which_seq in ['positive', 'negative']:
        fragments_od = load_fasta_gz_and_return_OrderedDict(\
                                '../../CNN_dir/100bp_seqs/100bp_train_test_fastas_dir/100bp_'+which_seq+'_test_chr'+chr_index+'.fasta.gz',\
                                input_len)
        scores_od = log_likelihood_each_matrix_global(PWM_fname, chr_index, fragments_od, input_len)

        with open('../loglikelihood_ratio_scores/'+which_seq+'_chr_'+chr_index+'_global_comp_made_motif'+PWM_num+'.json', 'w') as f:
            json.dump(scores_od,f)


if __name__ == '__main__':
    args = sys.argv
    input_len = 100

    for chr_index in [str(1),str(2),str(3),str(4),str(5),str(6),str(7),str(8),\
                      str(9),str(10),str(11),str(12),str(13),str(14),str(15),str(16), \
                      str(17),str(18),str(19),str(20),str(21),str(22),'X']:
    #for chr_index in ['N']:
        for PWM_num in [str(i) for i in range(16,20)]:
            main(chr_index=chr_index, PWM_num=PWM_num, input_len=input_len)