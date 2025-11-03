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
  ## ① 只选非灰色（有细胞类型颜色）的点
  filter(Celltype != "Multiple or Other") %>%      

  ## ② 仍然按显著性与效应筛选（阈值可自行调整）
  filter(FDR < 0.01, abs(logFC) > 1, !is.na(geneSymbol)) %>% 

  ## ③ FDR 最小的排在前面
  arrange(FDR) %>%                                  

  ## ④ 同一基因只保留一条
  distinct(geneSymbol, .keep_all = TRUE) %>%        

  ## ⑤ 最多 10 个标签
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
    type   = "open",                # 空心箭头
    length = unit(0.025, "npc"),    # ← 更小的箭头
    ends   = "last",
    angle  = 25
  )
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "log(FC)", y = "-log10(FDR)", color = "Cell Type") +
  theme_minimal()

