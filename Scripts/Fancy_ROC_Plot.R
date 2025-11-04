# ---- packages ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(pROC)
  library(tidyr)
  library(purrr)
})

# ---- parameters ----
set.seed(42)
top_n   <- 10        # how many genes to plot (by |AUC - 0.5|)
boot_n  <- 200       # bootstrap resamples for CI
grid_n  <- 100       # FPR grid points for ribbons
out_png <- "vitA_top_genes_ROC_R.png"
out_csv <- "vitA_top_genes_auc_R.csv"

# ---- load data ----
expr <- read_csv("data-raw/A-Retinol/Training/GSE299530_features_combined_2021.csv",
                 show_col_types = FALSE)

# Ensure the first column is named ENTREZID
names(expr)[1] <- "ENTREZID"
expr <- expr |>
  mutate(ENTREZID = suppressWarnings(as.integer(ENTREZID))) |>
  filter(!is.na(ENTREZID))

map_df <- read_csv("data-raw/mouse_entrez_to_symbol.cleaned.csv",
                   show_col_types = FALSE) |>
  transmute(ENTREZID = suppressWarnings(as.integer(ENTREZID)),
            SYMBOL   = as.character(SYMBOL)) |>
  filter(!is.na(ENTREZID)) |>
  distinct(ENTREZID, .keep_all = TRUE)

# ---- derive labels from column names ----
sample_cols <- setdiff(names(expr), "ENTREZID")

lab_vec <- case_when(
  str_detect(tolower(sample_cols), "vad") ~ 1L,   # positive class
  str_detect(tolower(sample_cols), "vas") ~ 0L,
  TRUE ~ NA_integer_
)

label_df <- tibble(sample = sample_cols, y = lab_vec) |>
  filter(!is.na(y))

labeled_cols <- label_df$sample
y            <- label_df$y

# Keep only labeled samples
X <- expr |>
  select(ENTREZID, all_of(labeled_cols))

# Convert to numeric and impute per-gene median for NAs
X[labeled_cols] <- lapply(X[labeled_cols], function(v) suppressWarnings(as.numeric(v)))
X <- X |>
  rowwise() |>
  mutate(across(all_of(labeled_cols), ~ ifelse(is.na(.x), median(c_across(all_of(labeled_cols)), na.rm = TRUE), .x))) |>
  ungroup()

# ---- compute per-gene AUC ----
calc_auc <- function(scores, y) {
  # if scores are constant, AUC is undefined wrt ranking
  if (isTRUE(all(scores == scores[1]))) return(NA_real_)
  as.numeric(pROC::auc(response = y, predictor = scores, quiet = TRUE))
}

auc_tbl <- X |>
  mutate(AUC = pmap_dbl(across(all_of(labeled_cols)), ~ calc_auc(c(...), y))) |>
  filter(!is.na(AUC)) |>
  mutate(Delta = abs(AUC - 0.5)) |>
  arrange(desc(Delta)) |>
  left_join(map_df, by = "ENTREZID") |>
  mutate(Label = if_else(!is.na(SYMBOL) & SYMBOL != "",
                         paste0(SYMBOL, " (AUC = ", sprintf("%.2f", AUC), ")"),
                         paste0(ENTREZID, " (AUC = ", sprintf("%.2f", AUC), ")")))

# Save a CSV of the top genes by |AUC-0.5|
auc_tbl |>
  slice_head(n = top_n) |>
  select(ENTREZID, SYMBOL, AUC) |>
  write_csv(out_csv)

# ---- helper: ROC with bootstrap CI ribbon ----
roc_with_ci <- function(scores, y, boot_n = 200, grid_n = 100) {
  # main ROC
  roc_obj <- pROC::roc(response = y, predictor = scores, quiet = TRUE)
  fpr     <- 1 - roc_obj$specificities
  tpr     <- roc_obj$sensitivities
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # common FPR grid to build ribbons
  grid <- seq(0, 1, length.out = grid_n)
  
  # bootstrap resampling
  n <- length(y)
  tpr_mat <- matrix(NA_real_, nrow = grid_n, ncol = boot_n)
  
  for (b in seq_len(boot_n)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    yb  <- y[idx]
    sb  <- scores[idx]
    # need at least one sample from each class
    if (length(unique(yb)) < 2) next
    roc_b <- pROC::roc(response = yb, predictor = sb, quiet = TRUE)
    fpr_b <- 1 - roc_b$specificities
    tpr_b <- roc_b$sensitivities
    # interpolate TPR at the common FPR grid
    tpr_interp <- approx(x = fpr_b, y = tpr_b, xout = grid, ties = mean, rule = 2)$y
    tpr_interp[1]  <- 0
    tpr_interp[grid_n] <- 1
    tpr_mat[, b] <- tpr_interp
  }
  
  ci_lower <- apply(tpr_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  ci_upper <- apply(tpr_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  list(
    curve = tibble(FPR = fpr, TPR = tpr),
    ribbon = tibble(FPR = grid, TPR_low = ci_lower, TPR_high = ci_upper),
    auc = auc_val
  )
}

# ---- build plotting data for the top N genes ----
top_genes <- auc_tbl |>
  slice_head(n = top_n)

roc_list <- pmap(
  list(scores_row = split(X[labeled_cols], X$ENTREZID)[as.character(top_genes$ENTREZID)],
       gid = top_genes$ENTREZID,
       label = top_genes$Label),
  function(scores_row, gid, label) {
    scores <- as.numeric(scores_row)
    rc     <- roc_with_ci(scores, y, boot_n = boot_n, grid_n = grid_n)
    list(
      gene_id = gid,
      label   = label,
      curve   = rc$curve |> mutate(Gene = label),
      ribbon  = rc$ribbon |> mutate(Gene = label),
      auc     = rc$auc
    )
  }
)

curve_df  <- bind_rows(lapply(roc_list, `[[`, "curve"))
ribbon_df <- bind_rows(lapply(roc_list, `[[`, "ribbon"))

# ---- plot ----
p <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", alpha = 0.6) +
  geom_ribbon(data = ribbon_df,
              aes(x = FPR, ymin = TPR_low, ymax = TPR_high, fill = Gene),
              alpha = 0.15, color = NA) +
  geom_line(data = curve_df,
            aes(x = FPR, y = TPR, color = Gene), linewidth = 1) +
  scale_x_continuous("False Positive Rate", limits = c(0, 1), expand = c(0,0)) +
  scale_y_continuous("True Positive Rate",  limits = c(0, 1), expand = c(0,0)) +
  guides(fill = "none") +
  labs(title = "ROC Curves for Vitamin A Signatures",
       subtitle = paste0("Shaded 95% bootstrap CI (n = ", boot_n, "); top ", top_n, " genes by |AUC - 0.5|"),
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(out_png, p, width = 10, height = 6, dpi = 160)

message("Saved figure: ", out_png)
message("Saved AUC table: ", out_csv)
