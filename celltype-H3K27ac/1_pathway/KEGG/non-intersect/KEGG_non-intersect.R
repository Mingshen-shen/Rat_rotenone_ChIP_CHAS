# Note: -log10(FDR) values may differ slightly (e.g., 0.01-0.02)
# due to changes in annotation databases or package versions.
# These differences do not affect overall interpretation.

# Load packages
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# Load TxDb
txdb <- TxDb.Rnorvegicus.UCSC.rn7.refGene

# Cell type & BED file paths
cell_types <- c("NeuN", "Nurr", "Olig", "Neg")
bed_files <- c("~/Desktop/NeuN_cleaned.bed",
               "~/Desktop/Nurr_cleaned.bed",
               "~/Desktop/Olig_cleaned.bed",
               "~/Desktop/Neg_cleaned.bed")

# Color palette (same as GO non-intersect)
cell_colors <- c(
  "Neg"  = "#E15759",
  "NeuN" = "#59A14F",
  "Nurr" = "#4E79A7",
  "Olig" = "#B07AA1"
)

# Data frame to store all results
all_kegg <- data.frame()

# KEGG analysis for each cell type
for (i in seq_along(cell_types)) {

  cell <- cell_types[i]
  bed <- bed_files[i]

  message("Processing ", cell, "...")

  # Peak annotation
  peak_annot <- annotatePeak(bed, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Rn.eg.db")

  # Extract gene IDs
  genes <- as.data.frame(peak_annot)$geneId

  # KEGG enrichment
  kegg_res <- enrichKEGG(gene = genes,
                         organism = 'rno',
                         pvalueCutoff = 0.05)

  # Convert to data frame and filter
  kegg_df <- as.data.frame(kegg_res) %>%
    filter(p.adjust <= 0.05) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)

  # Add cell type
  kegg_df$CellType <- cell

  # Combine results
  all_kegg <- bind_rows(all_kegg, kegg_df)
}

# Add label column (rounded to 2 decimal places)
all_kegg$label <- round(-log10(all_kegg$p.adjust), 2)

# Create dotplot with labels
ggplot(all_kegg, aes(x = CellType, y = reorder(Description, p.adjust))) +
  geom_point(aes(size = -log10(p.adjust), color = CellType)) +
  geom_text(aes(label = label), hjust = -0.4, size = 3) +
  scale_color_manual(values = cell_colors) +
  scale_size_continuous(name = "-log10(FDR)") +
  labs(x = "Cell Type", y = "KEGG Pathway",
       title = "KEGG Pathway Enrichment (Top 10, FDR <= 0.05)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
