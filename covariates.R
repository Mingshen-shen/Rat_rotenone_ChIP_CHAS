# Function: Using 4 cell types as covariates
group_vector_cortex <- pheno$Group  # or factor(pheno$Group) grouping vector


da_neuron_scores <- PD_MF_noBAM_SN$proportions[, c("Dopaminergic", "Neuron")] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Sample")

# Assuming your score table is named score_data
oligo_score_data <- score_data %>%
  filter(Celltype == "Oligodendrocyte") %>%
  select(Sample, Score) %>%
  rename(oligo_score = Score)

# Merge with pheno by Sample, ensuring consistent order.
pheno_with_oligo <- pheno %>%
  left_join(oligo_score_data, by = "Sample")

# Generate Covariate Vector
oligo_score_vector <- pheno_with_oligo$oligo_score

chas_cortex_scores <- celltype_scores_cortex %>%
  filter(Celltype %in% c("Dopaminergic", "Neg", "Neuron", "Oligodendrocyte")) %>%
  pivot_wider(names_from = Celltype, values_from = Score)

pheno_with_chas <- pheno %>%
  left_join(chas_cortex_scores, by = "Sample")
