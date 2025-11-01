group_vector <- factor(PD_pheno_SN$Group)

# Constructing a DGEList object
dge <- DGEList(counts = counts_mat)

# Filter out low-expression peaks (optional)
keep <- filterByExpr(dge, group = group_vector)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Standardization
dge <- calcNormFactors(dge)

# Constructing the design matrix
design <- model.matrix(~ group_vector)


# Estimated Dispersion
dge <- estimateDisp(dge, design)

# Fitting Model
fit <- glmFit(dge, design)

# Perform statistical testing: Test the effect of group_vector (column 2)
lrt <- glmLRT(fit, coef = 2)

# Extract significantly different peaks (e.g., FDR < 0.05)
diff_peaks <- topTags(lrt, n = Inf)
diff_peaks_df <- as.data.frame(diff_peaks)
sig_peaks <- diff_peaks_df[diff_peaks_df$FDR < 0.05, ]

# Further filtering can be performed based on logFC
hyper_peaks <- rownames(sig_peaks)[sig_peaks$logFC > 0]
hypo_peaks  <- rownames(sig_peaks)[sig_peaks$logFC < 0]
