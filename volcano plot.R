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
  ## ① （Optional）Select only non-gray points (those with cell type colors).
  filter(Celltype != "Multiple") %>%      

  ## ② Still filter based on significance and effect size (thresholds can be adjusted as needed)
  filter(FDR < 0.01, abs(logFC) > 1, !is.na(geneSymbol)) %>% 

  ## ③ FDR: Smallest first
  arrange(FDR) %>%                                  

  ## ④ Only one copy of the same gene is retained.
  distinct(geneSymbol, .keep_all = TRUE) %>%        

  ## ⑤ Up to 10 tags
  head(10)

# Volcano plot: Color-code cell types and highlight significant gene names.
ggplot(volcano_df, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = Celltype), alpha = .75, size = 2) +
  scale_color_manual(values = celltype_colors) +
  
 geom_text_repel(
  data = label_df,
  aes(label = geneSymbol),
  size = 4,
  max.overlaps = 15,
  segment.color = "black",
  segment.size  = .6,
  segment.alpha = .9,
  min.segment.length = 0,
  arrow = arrow(
    type   = "open",                
    length = unit(0.025, "npc"),    
    ends   = "last",
    angle  = 25
  )
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "log(FC)", y = "-log10(FDR)", color = "Cell Type") +
  theme_minimal()


