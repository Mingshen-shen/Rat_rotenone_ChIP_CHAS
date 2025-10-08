# Drawing volcano plot
dar_results$peakID <- rownames(dar_results)

volcano_df <- dar_results %>%
  left_join(peak_annot_SN, by = c("peakID" = "bulkPeak")) %>%  # add Celltype
  left_join(peak_annotation_df, by = c("peakID" = "peakID"))   # add gene names
# Process cell type annotation colors (classified into 4 specific types + Multiple/Other)
volcano_df$Celltype[is.na(volcano_df$Celltype)] <- "Multiple"

# Set color
celltype_colors <- c(
  "Neuron" = "#6BAED6",
  "Oligodendrocyte" = "#FBB4AE",
  "Dopaminergic" = "#FDAE6B",
  "Neg" = "#FC9272",
  "Multiple" = "gray70"
)
# Select significantly enriched genes (e.g., FDR < 0.01 and abs(logFC) > 1)
label_df <- volcano_df %>%
  
  filter(Celltype != "Multiple or Other") %>%      
  
  filter(FDR < 0.01, abs(logFC) > 1, !is.na(geneSymbol)) %>% 
  
  arrange(FDR) %>%                                  
  
  distinct(geneSymbol, .keep_all = TRUE) 



# Volcano plot: Color-code cell types and highlight significant gene names.
ggplot(volcano_df, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = Celltype), alpha = .75, size = 2) +
  scale_color_manual(values = celltype_colors) +
  
  geom_text_repel(
    data    = label_df,
    aes(label = geneSymbol),
    size    = 4,
    max.overlaps = 15,
    segment.color = "black",
    segment.size  = .6,
    segment.alpha = .9,
    min.segment.length = 0,
    arrow = arrow(
      type   = "open",           
      length = unit(0.02, "npc"),
      ends   = "last",           
      angle  = 25                
    )
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "log(FC)", y = "-log10(FDR)", color = "Cell Type") +
  theme_minimal()
