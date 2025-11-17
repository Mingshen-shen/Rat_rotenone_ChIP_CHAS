# Load the required objects. Replace with your own file path.
celltype_scores_SN <- readRDS("celltype_scores_SN.rds")
counts_mat <- readRDS("counts_mat.rds")

library(dplyr)

# 1. Here, taking the Oligodendrocyte in CHAS_score as an example, only the Sample and Score are retained.
oligo_score_df <- celltype_scores_SN %>%
  filter(Celltype == "Oligodendrocyte") %>%
  select(Sample, Score)

# 2. Ensure the sample order matches the column names in the count matrix
#     (Critical! Otherwise the design will be misaligned)
oligo_score_df <- oligo_score_df %>%
  arrange(match(Sample, colnames(counts_mat)))

# 3. Generate the covariate vector (sorted by sample order)
oligo_chas_scores <- oligo_score_df$Score


