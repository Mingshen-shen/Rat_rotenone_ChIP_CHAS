PD_consensusPeaks_SN <- ConsensusPeaks(bulk_SN_df, counts_mat, 
                                    bed_ref_df, counts_ref_matrix)

refSamples <- data.frame(
  Sample = colnames(counts_ref_matrix),
  CellType = c("Neuron", "Neuron", "Neuron", "Oligodendrocyte", "Oligodendrocyte", "Dopaminergic", "Dopaminergic", "Neg","Neg""),
  stringsAsFactors = FALSE
)
refSamples <- refSamples[order(refSamples$CellType), ]

PD_MF_noBAM_SN <- CelltypeProportion(
  PD_consensusPeaks_SN$newBulkCPM,
  PD_consensusPeaks_SN$newRefCPM,
  PD_consensusPeaks_SN$consensusPeaks,
  refSamples,
  signature = NULL
)

plot_MF_props(PD_MF_noBAM_SN, sampleLabel=FALSE)

plot_MF_groups(PD_MF_noBAM_SN, PD_pheno_SN)

plot_correlation(PD_MF_noBAM_SN, celltype_scores_SN, PD_pheno_SN)



