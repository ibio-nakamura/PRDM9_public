import json
from scipy import stats


CNN_spearman_dict = json.load(open('./record/CNN_spearman_data.json'))
PWM_spearman_dict = json.load(open('./record/PWM_spearman_data.json'))


PWM_cold_score_list = PWM_spearman_dict['cold_score_list']
PWM_cold_rec_list = PWM_spearman_dict['cold_rec_list']
CNN_cold_score_list = CNN_spearman_dict['cold_score_list']
CNN_cold_rec_list = CNN_spearman_dict['cold_rec_list']
print(len(PWM_cold_score_list))
print(len(CNN_cold_score_list))

print(CNN_cold_score_list[:10])
print(CNN_cold_score_list[-10:])
with open('./spearmanr/CNN_spearmanr_in_ChIP-negative.txt','w') as f:
    f.write(str(stats.spearmanr(CNN_cold_score_list,CNN_cold_rec_list)))

print(PWM_cold_score_list[:10])
print(PWM_cold_score_list[-10:])
with open('./spearmanr/PWM_spearmanr_in_ChIP-negative.txt','w') as f:
    f.write(str(stats.spearmanr(PWM_cold_score_list,PWM_cold_rec_list)))
