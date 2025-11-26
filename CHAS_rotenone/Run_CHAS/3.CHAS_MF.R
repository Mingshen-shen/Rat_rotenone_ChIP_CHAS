# Script to run CHAS-MF deconvolution of bulk susbtantia nigra H3K27ac profiles
# using a sorted cell reference stored in the objects/ directory.
# CHAS available on github: https://github.com/Marzi-lab/CHAS.

# Load R packages
library(CHAS)
library(ggplot2)

# Load the required objects. Replace with your own file path.
bed_ref_df <- readRDS("objects/bed_ref_df.rds")
counts_ref_matrix <- readRDS("objects/counts_ref_matrix.rds")
bulk_SN_df <- readRDS("objects/bulk_SN_df.rds")
counts_mat <- readRDS("objects/counts_mat.rds")
PD_pheno_SN <- readRDS("objects/PD_pheno_SN.rds")
celltype_scores_SN <- readRDS("objects/celltype_scores_SN.rds")

# Create consensus peaks
PD_consensusPeaks_SN <- ConsensusPeaks(bulk_SN_df, counts_mat, 
                                    bed_ref_df, counts_ref_matrix)

refSamples <- data.frame(
  Sample = colnames(counts_ref_matrix),
  CellType = c("Neuron", "Neuron", "Neuron", "Oligodendrocyte", "Oligodendrocyte", "Dopaminergic", "Dopaminergic", "Neg","Neg"),
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

pdf("CHAS_Results/SN_proportions.pdf", width = 6, height = 4)
plot_MF_props(PD_MF_noBAM_SN, sampleLabel=FALSE)
dev.off()

p <- plot_MF_groups(PD_MF_noBAM_SN, PD_pheno_SN)
ggsave("CHAS_Results/SN_MF.pdf", plot = p, width = 6, height = 4, units = "in", device = cairo_pdf)

p <- plot_correlation(PD_MF_noBAM_SN, celltype_scores_SN, PD_pheno_SN)
ggsave("CHAS_Results/correlation_SN.pdf", plot = p, width = 6, height = 4, units = "in", device = cairo_pdf)


