# Hypergeometric Analysis
# Load the required objects. Replace with your own file path.
peak_annot_SN <- readRDS("peak_annot_SN.rds")
universe_peaks_SN <- readRDS("universe_peaks_SN.rds")

# Required R packages
library(pheatmap)

# You should already have these datasets:
  
  # dar_status_vector: Status (Hyper, Hypo, NS) for all peak IDs

  # peak_annot_SN: Cell type to which each peak belongs

  # universe_peaks_SN: all the peaks

# Constructing a hypergeometric test function
run_hypergeometric <- function(dar_list, celltype, peak_annot, universe_peaks) {
  # A: Number of peaks belonging to this cell type in DAR
  A <- sum(dar_list %in% peak_annot$bulkPeak[peak_annot$Celltype == celltype])
  
  # B: The number of peaks belonging to this cell type across the entire universe
  B <- sum(peak_annot$bulkPeak %in% universe_peaks & peak_annot$Celltype == celltype)
  
  # C: Number of peaks in the entire universe that do not belong to this cell type
  C <- length(universe_peaks) - B
  
  # D: Number of peaks in DAR (total successful matches)
  D <- length(dar_list)
  
  # p-value: The probability that at least A peaks are randomly sampled from the universe, where D peaks are selected, are from this cell type
  pval <- phyper(q = A - 1, m = B, n = C, k = D, lower.tail = FALSE)
  
  return(pval)
}

# Perform enrichment analysis across multiple cell types and directions
# Set cell types
celltypes <- c("Neuron", "Oligodendrocyte", "Dopaminergic", "Neg")

# Results table
enrichment_matrix <- matrix(NA, nrow = 2, ncol = length(celltypes))
rownames(enrichment_matrix) <- c("Hyper", "Hypo")
colnames(enrichment_matrix) <- celltypes

# Hyper
for (ct in celltypes) {
  enrichment_matrix["Hyper", ct] <- run_hypergeometric(hyper_peaks, ct, peak_annot_SN, universe_peaks_SN)
}

# Hypo
for (ct in celltypes) {
  enrichment_matrix["Hypo", ct] <- run_hypergeometric(hypo_peaks, ct, peak_annot_SN, universe_peaks_SN)
}

# Convert to -log10(p) and perform FDR adjustment
logFDR_matrix <- -log10(apply(enrichment_matrix, c(1,2), function(p) p.adjust(p, method = "fdr")))

# Significance (FDR < 0.05)
sig_matrix <- enrichment_matrix < 0.05

# Draw heatmap
pheatmap(logFDR_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(50),
         display_numbers = ifelse(sig_matrix, "*", ""),
         fontsize_number = 20,
         main = "Hypergeometric Enrichment (-log10 FDR)") 

