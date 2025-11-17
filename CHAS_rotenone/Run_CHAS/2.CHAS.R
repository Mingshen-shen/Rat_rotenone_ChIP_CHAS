# Load the required objects. Replace with your own file path.
neg_peaks <- readRDS("neg_peaks.rds")
da_peaks <- readRDS("da_peaks.rds")
neuron_peaks <- readRDS("neuron_peaks.rds")
oligodendrocyte_peaks <- readRDS("oligodendrocyte_peaks.rds")
bulk_SN_df <- readRDS("bulk_SN_df.rds")
counts_mat <- readRDS("counts_mat.rds")
PD_pheno_SN <- readRDS("PD_pheno_SN.rds")


celltype_list <- list(Dopaminergic = da_peaks, Neg = neg_peaks, Neuron = neuron_peaks, Oligodendrocyte = oligodendrocyte_peaks)

celltype_specific_peaks_SN <- CelltypeSpecificPeaks(bulk_SN_df, celltype_list, 0.5)

celltype_scores_SN <- CelltypeScore(counts_mat, celltype_specific_peaks_SN, method="mean") 

plot_celltype_annotations(celltype_specific_peaks_SN)


plot_celltype_scores(celltype_scores_SN, PD_pheno_SN)

