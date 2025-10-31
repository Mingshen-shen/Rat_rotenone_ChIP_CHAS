celltype_list <- list(Dopaminergic = da_peaks, Negative = neg_peaks, Neuron = neuron_peaks, Oligodendrocyte = oligodendrocyte_peaks)

celltype_specific_peaks_SN <- CelltypeSpecificPeaks(bulk_SN_df, celltype_list, 0.5)

celltype_scores_SN <- CelltypeScore(counts_mat, celltype_specific_peaks_SN, method="mean") 

plot_celltype_annotations(celltype_specific_peaks_SN)

plot_celltype_scores(celltype_scores_SN, PD_pheno_SN)