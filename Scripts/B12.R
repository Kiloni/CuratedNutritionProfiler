
# ============================
# FULL PIPELINE FOR B12 ANALYSIS IN R
# ============================

# Load required packages
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ROCR)
library(caret)
library(stats)
library(reshape2)

# ----------------------------
# STEP 1: Load and Merge Datasets
# ----------------------------
files <- c(
  "~/Downloads/GSE103417_series_matrix.txt",
  "~/Downloads/GSE104164_series_matrix.txt",
  "~/Downloads/GSE161763_series_matrix.txt",
  "~/Downloads/GSE99113_JC_NS190_AstrocyteB12.txt",
  "~/Downloads/GSE104164_MDD-rat-cerebellum_normalized-data_with_annotation.txt"
)

dfs <- lapply(files, function(f) {
  read.delim(f, header = TRUE, comment.char = "!")
})

# Merge on union of genes
merged_df <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE), dfs)
rownames(merged_df) <- merged_df$Row.names
merged_df <- merged_df[, -1]

cat("Merged dataset dimensions:", dim(merged_df), "\n")
cat("Missing values:", sum(is.na(merged_df)), "\n")

# ----------------------------
# STEP 2: Normalization Check
# ----------------------------

# Identify numeric columns
numeric_cols <- sapply(merged_df, is.numeric)

# Apply log2 normalization only to numeric columns
if (max(merged_df[, numeric_cols], na.rm = TRUE) > 100) {
  merged_df[, numeric_cols] <- log2(merged_df[, numeric_cols] + 1)
}

# ----------------------------
# STEP 3: PCA Visualization
# ----------------------------
# Transpose for PCA
pca_data <- t(merged_df)
pca_res <- prcomp(na.omit(pca_data), scale. = TRUE)

# Keep only numeric columns for PCA
numeric_cols <- sapply(merged_df, is.numeric)
pca_data <- t(merged_df[, numeric_cols])

# Run PCA
pca_res <- prcomp

``

## 1) Convert matrix to data.frame
pca_df_in <- as.data.frame(pca_data)

## 2) Replace Inf/-Inf with NA to treat uniformly
pca_df_in[!is.finite(as.matrix(pca_df_in))] <- NA

## 3) Drop rows with any NA and columns with zero variance
row_ok <- apply(pca_df_in, 1, function(r) all(is.finite(r)))
pca_clean <- pca_df_in[row_ok, , drop = FALSE]

const_cols <- which(apply(pca_clean, 2, function(x) var(x, na.rm = TRUE) == 0))
if (length(const_cols) > 0) {
  pca_clean <- pca_clean[, -const_cols, drop = FALSE]
}

## 4) Run PCA (center + scale recommended when features have different magnitudes)
pca_res <- stats::prcomp(pca_clean, center = TRUE, scale. = TRUE)

## 5) Get scores robustly and align B12_status to filtered rows
scores <- predict(pca_res)
B12_status_clean <- B12_status[row_ok]
stopifnot(nrow(scores) == length(B12_status_clean))

## 6) Final frame for plotting
pca_df <- data.frame(
  PC1 = scores[, 1],
  PC2 = scores[, 2],
  B12_status = B12_status_clean
)

head(pca_df)


# Simulate B12 status (replace with actual annotation column)
B12_status <- factor(rep(c("Deficient", "Sufficient"), length.out = nrow(pca_data)))

pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], B12_status = B12_status)

ggplot(pca_df, aes(x = PC1, y = PC2, color = B12_status)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA by B12 Status")

# ----------------------------
# STEP 4: Heatmap of Top Variable Genes
# ----------------------------
gene_var <- apply(merged_df, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:50]

pheatmap(merged_df[top_genes, ], scale = "row", color = colorRampPalette(c("blue", "white", "red"))(50))

# ----------------------------
# STEP 5: Enrichment Analysis (Placeholder)
# ----------------------------
cat("Run enrichment using clusterProfiler or gprofiler2 here.\n")

# ----------------------------
# STEP 6: Gene Signature Analysis
# ----------------------------
# ANOVA for all genes
anova_results <- sapply(rownames(merged_df), function(gene) {
  df <- data.frame(expr = as.numeric(merged_df[gene, ]), status = B12_status)
  fit <- aov(expr ~ status, data = df)
  summary(fit)[[1]][["Pr(>F)"]][1]
})

anova_df <- data.frame(Gene = names(anova_results), ANOVA_p = anova_results)
anova_df <- anova_df[order(anova_df$ANOVA_p), ]

# ROC/AUC function
plot_roc_for_genes <- function(genes, title) {
  pred_df <- data.frame()
  for (gene in genes) {
    y_true <- ifelse(B12_status == "Sufficient", 1, 0)
    y_score <- as.numeric(merged_df[gene, ])
    pred <- prediction(y_score, y_true)
    perf <- performance(pred, "tpr", "fpr")
    auc_val <- performance(pred, "auc")@y.values[[1]]
    lines(perf@x.values[[1]], perf@y.values[[1]], col = sample(colors(), 1), lwd = 2)
    pred_df <- rbind(pred_df, data.frame(Gene = gene, AUC = auc_val))
  }
  abline(0, 1, lty = 2)
  title(main = title)
  legend("bottomright", legend = paste(pred_df$Gene, "AUC=", round(pred_df$AUC, 2)), cex = 0.8)
}

# Plot ROC for top genes
top_genes_anova <- head(anova_df$Gene, 5)
plot_roc_for_genes(top_genes_anova, "Top Predictive Genes ROC")

# Provided gene list
B12_Cobalamin <- c("E2F1", "Rps3", "Rps6", "Rpl5", "Rpl7", "Rpl11", "Rpl13a", "Rps14", "Rps27a", "Rpl22")
plot_roc_for_genes(B12_Cobalamin, "B12_Cobalamin Gene List ROC")
``
