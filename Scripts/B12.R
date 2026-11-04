# =============================================================================
# B12 ANALYSIS — rewritten 2026-11-04
#
# What changed from the previous version (still in git history):
#   * It read four files from ~/Downloads that no longer exist. It now reads the
#     processed tables in data/processed/Vitamin_B12/ (made by
#     Scripts/process_cnp_source_data.R from fresh GEO downloads).
#   * It used SIMULATED labels: B12_status <- rep(c("Deficient","Sufficient"), ...)
#     alternated sample by sample, so its ANOVA, ROC curves and any gene lists or
#     Enrichr results built from them do not reflect B12 status. Labels now come
#     from each study's GEO annotation (the *_pheno.csv files).
#   * It merged mouse arrays, mouse RNA-seq and rat arrays by row name and analysed
#     them as one matrix. Each dataset is now analysed on its own, since platforms and
#     species are not directly comparable.
#
# Per dataset: PCA coloured by B12 status, a heatmap of the B12_Cobalamin genes,
# per-gene ANOVA (adjusted for the study covariate where there is one) and per-gene
# ROC AUC. Signature-level validation (GSVA/ssGSEA/PLAGE/zScore, LR, RF, LOOCV and
# external testing) is in Reports/CNP_Signature_Validation.Rmd.
#
# Outputs: Reports/results/B12/  (B12_gene_statistics.csv, B12_plots.pdf)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(pheatmap); library(pROC)
})
PROJECT_DIR <- if (file.exists("curatedNutritionDataset.Rproj")) getwd() else
  path.expand("~/RStudio Projects/curatedNutritionDataset")
OUT <- file.path(PROJECT_DIR, "Reports", "results", "B12"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ex <- parse(file.path(PROJECT_DIR, "Signatures", "NutSignaturesList.R")); eval(ex[[1]])
b12_genes <- NutSignatures$B12_Cobalamin
rodent_case <- function(g) paste0(substr(g, 1, 1), tolower(substring(g, 2)))   # E2F1 -> E2f1

# Case = B12-deficient / B12-deprived; control = replete (same definitions as the report)
datasets <- list(
  list(id = "GSE161763", case = "cblX_F80L",          control = "Control",                    covar = NULL,     label = "Mouse brain, cblX (Ronin F80L)"),
  list(id = "GSE103417", case = "TO_B12_sequestered", control = c("WT", "OT_control_fusion"),  covar = "medium", label = "Mouse N1E-115, B12 sequestration"),
  list(id = "GSE104164", case = "MDD",                control = "Control",                    covar = NULL,     label = "Rat cerebellum, methyl-donor deficiency"),
  list(id = "GSE99113",  case = "B12_deficient",      control = "B12_sufficient",             covar = NULL,     label = "Mouse astrocytes, B12-free medium (2 vs 2)")
)

stats <- list()
pdf(file.path(OUT, "B12_plots.pdf"), width = 9, height = 6)
for (d in datasets) {
  e <- fread(file.path(PROJECT_DIR, "data", "processed", "Vitamin_B12", paste0(d$id, "_expr.csv.gz")), data.table = FALSE)
  m <- as.matrix(e[, -1]); rownames(m) <- e$gene
  p <- read.csv(file.path(PROJECT_DIR, "data", "processed", "Vitamin_B12", paste0(d$id, "_pheno.csv")))
  p$status <- ifelse(p$condition %in% d$case, "Deficient", ifelse(p$condition %in% d$control, "Sufficient", NA))
  p <- p[!is.na(p$status), ]; m <- m[, p$sample]

  # PCA on the 2,000 most variable genes
  top <- head(order(apply(m, 1, var), decreasing = TRUE), 2000)
  pc <- prcomp(t(m[top, ]), scale. = TRUE); ve <- round(100 * summary(pc)$importance[2, 1:2])
  print(ggplot(data.frame(pc$x[, 1:2], p), aes(PC1, PC2, colour = status, shape = condition)) +
          geom_point(size = 3) + theme_minimal() +
          labs(title = paste(d$id, "-", d$label), x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")))

  # B12_Cobalamin genes present in this dataset
  g <- intersect(rodent_case(b12_genes), rownames(m))
  if (length(g) >= 2) pheatmap(m[g, ], scale = "row", annotation_col = data.frame(status = p$status, row.names = p$sample),
                               main = paste(d$id, "- B12_Cobalamin genes (row z-score)"), fontsize_row = 8)

  # Per-gene ANOVA (covariate entered first) and ROC AUC for deficiency
  cv <- d$covar[vapply(d$covar, function(v) length(unique(p[[v]])) > 1, TRUE)]
  for (gene in g) {
    df <- data.frame(expr = m[gene, ], status = factor(p$status), p)
    f <- as.formula(paste("expr ~", paste(c(cv, "status"), collapse = " + ")))
    a <- summary(aov(f, data = df))[[1]]; r <- a[trimws(rownames(a)) == "status", ]
    auc <- as.numeric(auc(roc(df$status == "Deficient", df$expr, levels = c(FALSE, TRUE), direction = "<", quiet = TRUE)))
    stats[[length(stats) + 1]] <- data.frame(dataset = d$id, label = d$label, gene = gene,
      mean_deficient = mean(df$expr[df$status == "Deficient"]), mean_sufficient = mean(df$expr[df$status == "Sufficient"]),
      F = r[["F value"]], p_value = r[["Pr(>F)"]], AUC_higher_in_deficient = auc,
      n_deficient = sum(df$status == "Deficient"), n_sufficient = sum(df$status == "Sufficient"))
  }
}
dev.off()

res <- bind_rows(stats) |> group_by(dataset) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup()
fwrite(res, file.path(OUT, "B12_gene_statistics.csv"))
print(res |> mutate(across(where(is.numeric), ~ signif(.x, 3))) |> as.data.frame())
message("Saved: ", file.path(OUT, "B12_gene_statistics.csv"), " and B12_plots.pdf")
