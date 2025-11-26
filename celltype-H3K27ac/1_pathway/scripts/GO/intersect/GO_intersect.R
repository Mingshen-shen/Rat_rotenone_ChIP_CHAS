# BED file paths for each cell type
bed_files <- list(
    NeuN = "~/Desktop/unique_peaks/NeuN_unique.bed",
    Nurr = "~/Desktop/unique_peaks/Nurr_unique.bed",
    Olig = "~/Desktop/unique_peaks/Olig_unique.bed",
    Neg  = "~/Desktop/unique_peaks/Neg_unique.bed"
)

# GO enrichment analysis function
run_go <- function(bed_path, cell_type) {
    peak_annot <- annotatePeak(bed_path, TxDb = txdb, annoDb = "org.Rn.eg.db")
    genes <- na.omit(as.data.frame(peak_annot)$geneId)

    ego_bp <- enrichGO(gene = genes, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "BP", qvalueCutoff = 0.05) %>%
        as.data.frame() %>% mutate(Category = "BP", CellType = cell_type)

    ego_cc <- enrichGO(gene = genes, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "CC", qvalueCutoff = 0.05) %>%
        as.data.frame() %>% mutate(Category = "CC", CellType = cell_type)

    ego_mf <- enrichGO(gene = genes, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "MF", qvalueCutoff = 0.05) %>%
        as.data.frame() %>% mutate(Category = "MF", CellType = cell_type)

    bind_rows(ego_bp, ego_cc, ego_mf)
}

# Run analysis for each cell type
go_all <- bind_rows(
    run_go(bed_files$NeuN, "NeuN"),
    run_go(bed_files$Nurr, "Nurr"),
    run_go(bed_files$Olig, "Olig"),
    run_go(bed_files$Neg,  "Neg")
)

# Filter significant GO terms and select Top 10
go_sig <- go_all %>% filter(p.adjust < 0.05)
go_top <- go_sig %>% group_by(Category) %>% slice_min(order_by = p.adjust, n = 10) %>% ungroup()

# Plot
ggplot(go_top, aes(x = -log10(p.adjust), y = Description,
                   shape = Category, color = CellType)) +
    geom_point(size = 5) +
    labs(title = "Top 10 GO Terms per Category",
         x = "-log10(FDR)", y = "GO Term") +
    theme_minimal()
