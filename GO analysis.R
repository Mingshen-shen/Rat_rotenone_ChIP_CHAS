# Filter out hypo peaks that appear in the annotation.
hypo_annot <- peak_annot_SN[peak_annot_SN$bulkPeak %in% hypo_peaks, ]

# Retain only the hypo annotations, group them by cell type, and extract the peak IDs.
celltype_specific_peaks_hypo <- split(hypo_annot$bulkPeak, hypo_annot$Celltype)

# Ensure that no duplicates exist within each vector.
celltype_specific_peaks_hypo <- lapply(celltype_specific_peaks_hypo, unique)

# You can create a GRanges object for annotation:
# Assuming you have the following dataframe
# peak_coords <- data.frame(chr = ..., start = ..., end = ..., peakID = ...)

gr_peaks <- GRanges(seqnames = peak_coords$chr,
                    ranges = IRanges(start = peak_coords$start,
                                     end = peak_coords$end),
                    peakID = peak_coords$peakID)

# annotating peaks
peak_annot_SN <- annotatePeak(gr_peaks,
                              TxDb = TxDb.Rnorvegicus.UCSC.rn6.refGene,
                              annoDb = "org.Rn.eg.db")

# Convert to a dataframe and extract the necessary information.
peak_annot_df <- as.data.frame(peak_annot)
peak_annot_df <- peak_annot_df[, c("peakID", "geneSymbol")]

#Build a list mapping cell type to gene symbol
hypo_gene_list <- lapply(celltype_specific_peaks_hypo, function(peaks) {
  genes <- peak_annot_df %>%
    filter(peakID %in% peaks) %>%
    pull(geneSymbol) %>%
    unique()
  genes[!is.na(genes) & genes != ""]
})

# Convert gene symbol to ENTREZID
hypo_entrez_list <- lapply(hypo_gene_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)$ENTREZID
})

# GO enrichment analysis (only performing hypo here) simultaneously converted to a dataframe
go_hyper_list_SN <- lapply(names(hyper_entrez_list_SN), function(celltype) {
  genes <- hyper_entrez_list_SN[[celltype]]
  enrichGO(gene = genes,
           OrgDb = org.Rn.eg.db,
           keyType = "ENTREZID",
           ont = "ALL",
           pAdjustMethod = "fdr",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2) %>%
    as.data.frame() %>%
    mutate(celltype = celltype)
})

go_hyper_SN_merged <- do.call(rbind, go_hyper_list_SN)


# Add Bulk GO analysis

bulk_hypo_genes <- peak_annotation_df %>%
  filter(peakID %in% hypo_peaks_SN) %>%
  pull(geneSymbol) %>%
  unique()
bulk_hypo_genes <- bulk_hypo_genes[!is.na(bulk_hypo_genes) & bulk_hypo_genes != ""]

bulk_hypo_entrez <- bitr(bulk_hypo_genes,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Rn.eg.db)$ENTREZID

bulk_go <- enrichGO(gene = bulk_hypo_entrez,
                    OrgDb = org.Rn.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2) %>%
  mutate(celltype = "Bulk")

bulk_go_df <- as.data.frame(bulk_go_SN_hypo)
bulk_go_df$celltype <- "Bulk"

# Merging bulk and cell types

go_hyper_SN <- rbind(go_hyper_df_SN, bulk_go_df_SN_hyper)

# Select the top 5 most significant pathways for each cell type

go_top_hyper <- go_hyper_df_SN %>%           
  arrange(p.adjust) %>%                      
  group_by(celltype) %>%                     
  slice_head(n = 5) %>%                      
  ungroup()

# Include FDR-characterized size and GO category as graphical markers
go_top$GO_shape <- ifelse(go_top$ONTOLOGY == "BP", 16,
                          ifelse(go_top$ONTOLOGY == "CC", 17, 15))  


# Label
go_top_SN_hypo$GO_class <- recode(go_top_SN_hypo$ONTOLOGY,
                                  BP = "Biological Process",
                                  CC = "Cellular Component",
                                  MF = "Molecular Function")

go_top_SN_hypo$logFDR <- -log10(go_top_SN_hypo$p.adjust)

ggplot(go_top_SN_hypo, aes(x = celltype, y = Description)) +
  geom_point(aes(size = logFDR, shape = GO_class, color = celltype)) +
  scale_shape_manual(values = c("Biological Process" = 16,
                                "Cellular Component" = 17,
                                "Molecular Function" = 15)) +
  scale_size_continuous(range = c(2, 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) +
  labs(title = "GO enrichment of celltype-specific hypoacetylated peaks",
       x = "Cell Type",
       y = "GO Term",
       size = "FDR (-log10)",
       shape = "GO Category")

