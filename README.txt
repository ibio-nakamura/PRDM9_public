The structure of the repository is as follows.
The hg19 genome data (chr*.fa.gz) was downloaded from the UCSC genome browser (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/) and stored under make_Noverlap_fragments/hg19/ with the unzipped version (chr*.fa).

Supplementary_materials/
	├ PRDM9_sm.pdf: Hyperparameter setting and metrics of other CNN-training methods are shown.
	└ motifs.xlsx: PWMs extracted from the trained EBCN are included.

make_Noverlap_fragments/
	└ 100bp_NoverlapQuarry.py: Create a 100bp fragment from hg19 and save it in Fasta format.


calc_rec_rate_starts_from_0/
	└ calc_100bp_Noverlap_rec_rate_starts_from_0.py: Calculate the recombination rate of 100 bp.


overall_plots/
	├ N_bp_all_experiment_subplots.py: Illustrate the correlation between the recombination rate and the score.
	├ fast_N_bp_all_experiment_subplots.py: Quickly illustrate from data recorded by N_bp_all_experiment_subplots.py.
	├ corr_calc.py: Calculate the recombination rate and score Spearman correlation.
	└ corr_calc_ChIP-seq-negative.py: Calculate the Spearman correlation between predicted score and recombination rate to ChIP-seq-negative fragment.


CNN_dir/
├ train_original_for_100bp_with_optuna.py: Adjust hyperparameters using positive data with peaks centered, train, and draw ROCs using test data.
├ train_allocated_for_100bp_with_optuna.py: Adjust hyperparameters using positive data with randomly displaced peaks, train, and draw ROCs using test data.
├ train_3aug_for_100bp_with_optuna.py: Adjusting hyperparameters using augmented data, training, and drawing ROCs using test data.
├ SavedModel_predict_testdata.py: Predict test data using a trained model.
├ SavedModel_predict_Noverlap_100bp_original.py: Predict the 100bp fragment using the trained model in train_original_for_100bp_with_optuna.py.
├ SavedModel_predict_Noverlap_100bp_allocated.py: Predict the 100bp fragment using the trained model in train_allocated_for_100bp_with_optuna.py.
├ SavedModel_predict_Noverlap_100bp_3aug.py: Predict the 100bp fragment using the trained model in train_3aug_for_100bp_with_optuna.py.
├ get_optimal_threshold.py: Determine the threshold for discriminating positive and negative from ROC.
├ plot_distance_btw_peak_and_CNN_positive_random.py: Output histograms of the distances between CNN-positive and ChIP-seq peaks, and between CNN-negative and ChIP-seq peaks, and store the information.
└ local/100bp_seqs/
	└ make_train_data.py: Create training and test data.


PWM_dir/
├ global_comp_each_testdata_calc_at_26/
│	├ judge_optimal_PWM_num_with_AUC.py: Output AUC to determine how many PWMs should be used.
│	├ integrate_scores_and_output_jsons.py: Integrate the loglikelihood ratio score to the 100 bp test fragment calculated for each PWM.
│	├ PWM_AUROC.py: Draw the ROC and save the information.
│	└ calculate_loglikelihood/
│		└ automatically_scan_test_data.sh: Perform scoring on test data.
│
└ global_comp_each_Noverlap_calc_at_34/
	├ integrate_scores_and_output_jsons.py: Summarize the loglikelihood-ratio-score to the 100 bp fragment calculated for each PWM.
	├ get_optimal_threshold.py: Determine the threshold for discriminating positive and negative from ROC.
	├ plot_distance_btw_peak_and_PWM_positive_random.py: The distances between PWM-positive and ChIP-seq peaks and PWM-negative and ChIP-seq peaks are output as histograms to store the information.
	├ fast_hist_plot.py: Quickly plot a histogram of distances using stored information.
	└ calculate_loglikelihood/
		└ automatically_scan_each_Noverlap.sh: Perform scoring to 100 bp fragments using PWM.


plot_distance/
	├ fast_hist_plot.py: Output a histogram of the distance between the ChIP-seq peak and the fragment.
	└ test_distances.py: Statistically test the distance between the ChIP-seq peak and the fragment.


other_metrics_new/
├ other_metrics.py: Output the prediction accuracy of the model in various metrics.
├ kmer_count.py: Perform statistical tests of the GC content in the test data.
└ 13mer_judge: Perform a 13-mer motif match.


motif_finding/
├ conv_each_filter.py: Perform scoring on random DNA sequences using filters in the convolutional layer of a trained CNN.
├ make_PWM.py: CNN: Create PWMs based on convolutional scores.
├ AUC_for_each_motif.py: The AUC is output from the loglikelihood-score to the test data by each PWM.
├ logo_make.py: Create a sequence logo from PWM.
├ global_comp_each_motifs_from_26
│	└automatically_scan_test_data.sh: Score the test data in each PWM.
│
└ made_PWMs
	└ dev_multi.py: Split the created multi_PWM.txt (use this before using AUC_for_each_motif.py).
