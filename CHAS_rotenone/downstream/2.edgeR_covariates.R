# Load the required objects. Replace with your own file path
PD_pheno_SN <- readRDS("PD_pheno_SN.rds")

library(edgeR)

group_vector <- factor(PD_pheno_SN$Group)

# Constructing the design matrix
design <- model.matrix(~ group_vector + oligo_chas_scores)

# Constructing a DGEList object
dge <- DGEList(counts = counts_mat)

# Standardization
dge <- calcNormFactors(dge)

# Filter out low-expression peaks (optional)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Estimated Dispersion
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)  # coef = 2 → GroupRotenone vs PBS

# Extract significantly different peaks (e.g., FDR < 0.05)
dar_results <- topTags(qlf, n = Inf)$table

# Adding label
dar_results$Status <- ifelse(dar_results$FDR < 0.05 & dar_results$logFC > 0, "Hyper",
                       ifelse(dar_results$FDR < 0.05 & dar_results$logFC < 0, "Hypo", "NS"))


# Further filtering can be performed based on logFC
hyper_peaks <- rownames(dar_results[dar_results$Status == "Hyper", ])
hypo_peaks  <- rownames(dar_results[dar_results$Status == "Hypo", ])

