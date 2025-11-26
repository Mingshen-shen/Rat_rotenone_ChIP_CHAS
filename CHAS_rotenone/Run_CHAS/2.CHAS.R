# Script to run CHAS deconvolution of bulk susbtantia nigra H3K27ac profiles
# using a sorted cell reference stored in the objects/ directory.
# CHAS available on github: https://github.com/Marzi-lab/CHAS.

# Load R packages
library(CHAS)
library(ggplot2)

# Load the required objects. Replace with your own file path.
neg_peaks <- readRDS("objects/neg_peaks.rds")
da_peaks <- readRDS("objects/da_peaks.rds")
neuron_peaks <- readRDS("objects/neuron_peaks.rds")
oligodendrocyte_peaks <- readRDS("objects/oligodendrocyte_peaks.rds")
bulk_SN_df <- readRDS("objects/bulk_SN_df.rds")
counts_mat <- readRDS("objects/counts_mat.rds")
PD_pheno_SN <- readRDS("objects/PD_pheno_SN.rds")

# Create reference cell type list
celltype_list <- list(Dopaminergic = da_peaks, Neg = neg_peaks, Neuron = neuron_peaks, Oligodendrocyte = oligodendrocyte_peaks)
# Find cell type specific peaks in bulk dataset  
celltype_specific_peaks_SN <- CelltypeSpecificPeaks(bulk_SN_df, celltype_list, 0.5)
# Calculate cell type scores
celltype_scores_SN <- CelltypeScore(counts_mat, celltype_specific_peaks_SN, method="mean") 

# Plot cell type annotations
p <- plot_celltype_annotations(celltype_specific_peaks_SN)
# Save plot
ggsave("CHAS_Results/SN_annotation.pdf", plot = p, width = 6, height = 4, units = "in", device = cairo_pdf)


# Plot cell type scores
p <- plot_celltype_scores(celltype_scores_SN, PD_pheno_SN)
# Save plot
ggsave("CHAS_Results/SN_betweengroups.pdf", plot = p, width = 6, height = 4, units = "in", device = cairo_pdf)

