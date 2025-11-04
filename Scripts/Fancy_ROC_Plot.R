## If you do not already have a gene signature, you can use a dataset to find the signature using this code:

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



-----



----

## If you DO already have a gene signature, then try this instead:

  # ---- Packages ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(pROC)    # for ROC/AUC
  library(tidyr)
})

# ---- Inputs ----
expr_file <- "data-raw/A-Retinol/Training/GSE299530_features_combined_2021.csv"
map_file  <- "data-raw/mouse_entrez_to_symbol.cleaned.csv"
out_png   <- "Figures/vitA_selected_genes_ROC.png"
out_csv   <- "Signatures/vitA_selected_genes_auc.csv"

# --- List your genes here ---
genes_requested <- c(
  "Lrat","Tnni3k","Dancr","Suv39h2","Apol7c","Tnni3",
  "Fndc5","St6galnac5","Neil3","Cxcr6","Slc47a1",
  "Myoz2","Cenpi","Slc17a7","Lrch2","Shcbp1l",
  "Cox6a2","Hspb7","Grik4","Gm7511"
)

# ---- Load data ----
map_df <- read_csv(map_file, show_col_types = FALSE) |>
  mutate(ENTREZID = suppressWarnings(as.integer(ENTREZID))) |>
  filter(!is.na(ENTREZID)) |>
  distinct(SYMBOL, .keep_all = TRUE)

expr0 <- read_csv(expr_file, show_col_types = FALSE)
colnames(expr0)[1] <- "ENTREZID"
expr0 <- expr0 |>
  mutate(ENTREZID = suppressWarnings(as.integer(ENTREZID))) |>
  filter(!is.na(ENTREZID))

# Remove sample columns that are all NA
all_na_cols <- names(expr0)[-1][sapply(expr0[-1], function(x) all(is.na(x)))]
expr <- expr0 |>
  select(-any_of(all_na_cols))

# ---- Labels from sample names ----
sample_cols <- setdiff(names(expr), "ENTREZID")
lab <- tibble(sample = sample_cols) |>
  mutate(label = case_when(
    str_detect(tolower(sample), "vad") ~ 1L,
    str_detect(tolower(sample), "vas") ~ 0L,
    TRUE ~ NA_integer_
  )) |>
  filter(!is.na(label))

stopifnot(nrow(lab) > 0)
y <- lab$label

# ---- Map SYMBOL -> ENTREZ and keep those present in expr ----
sym2ent <- map_df |>
  filter(SYMBOL %in% genes_requested) |>
  select(SYMBOL, ENTREZID)

present_entrez <- intersect(sym2ent$ENTREZID, expr$ENTREZID)
sel_map <- sym2ent |>
  filter(ENTREZID %in% present_entrez)

if (nrow(sel_map) == 0) stop("None of the requested genes were found in the expression matrix.")

# ---- Build gene-by-sample matrix ----
X <- expr |>
  select(ENTREZID, all_of(lab$sample)) |>
  arrange(ENTREZID)

# median-impute NAs per gene
X[lab$sample] <- apply(X[lab$sample], 1, function(row) {
  med <- suppressWarnings(median(row, na.rm = TRUE))
  if (is.infinite(med)) med <- 0
  ifelse(is.na(row), med, row)
}) |> t() |> as.data.frame()
X <- bind_cols(expr["ENTREZID"], X)
names(X) <- c("ENTREZID", lab$sample)

# ---- Helper: ROC with bootstrap TPR bands ----
interp_vec <- function(x, y, xout) {
  # linear interpolation like numpy.interp
  approx(x, y, xout, ties = "ordered", rule = 2)$y
}

compute_gene_roc <- function(sym, ent, scores_mat, y, fpr_grid = seq(0, 1, length.out = 75), B = 200L, seed = 42L) {
  s <- as.numeric(scores_mat[scores_mat$ENTREZID == ent, lab$sample])
  if (length(unique(s)) < 2) return(NULL)
  
  # pROC expects response first, then predictor
  roc_obj <- roc(response = y, predictor = s, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(auc(roc_obj))
  
  # Obtain raw ROC points
  # coords with thresholds give specificity; convert to FPR and TPR
  coords_df <- coords(roc_obj, ret = c("specificity","sensitivity"), transpose = FALSE)
  fpr <- 1 - coords_df$specificity
  tpr <- coords_df$sensitivity
  
  # Bootstrap TPR bands over a fixed FPR grid
  set.seed(seed)
  n <- length(y)
  tprs <- matrix(NA_real_, nrow = 0, ncol = length(fpr_grid))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    yb <- y[idx]; sb <- s[idx]
    if (length(unique(yb)) < 2) next
    roc_b <- roc(response = yb, predictor = sb, quiet = TRUE, direction = "<")
    cd <- coords(roc_b, ret = c("specificity","sensitivity"), transpose = FALSE)
    fpr_b <- 1 - cd$specificity
    tpr_b <- cd$sensitivity
    tpr_i <- interp_vec(fpr_b, tpr_b, fpr_grid)
    tpr_i[1] <- 0
    tprs <- rbind(tprs, tpr_i)
  }
  lo <- apply(tprs, 2, function(z) quantile(z, 0.025, na.rm = TRUE))
  hi <- apply(tprs, 2, function(z) quantile(z, 0.975, na.rm = TRUE))
  
  list(
    sym = sym,
    ent = ent,
    auc = auc_val,
    fpr = fpr,
    tpr = tpr,
    grid = fpr_grid,
    lo = lo,
    hi = hi
  )
}

fpr_grid <- seq(0, 1, length.out = 75)
res_list <- pmap(
  list(sel_map$SYMBOL, sel_map$ENTREZID),
  ~ compute_gene_roc(..1, ..2, X, y, fpr_grid = fpr_grid, B = 200L, seed = 42L)
)

res_list <- compact(res_list)

# ---- Collect for plotting ----
curve_df <- map_dfr(res_list, function(r) {
  tibble(SYMBOL = r$sym, FPR = r$fpr, TPR = r$tpr)
})

band_df <- map_dfr(res_list, function(r) {
  tibble(SYMBOL = r$sym, FPR = r$grid, LO = r$lo, HI = r$hi)
})

auc_df <- map_dfr(res_list, function(r) tibble(SYMBOL = r$sym, ENTREZID = r$ent, AUC = r$auc)) |>
  arrange(SYMBOL)

# ---- Plot ----
p <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50", alpha = 0.6) +
  geom_ribbon(data = band_df, aes(x = FPR, ymin = LO, ymax = HI, fill = SYMBOL), alpha = 0.15, color = NA) +
  geom_path(data = curve_df, aes(x = FPR, y = TPR, color = SYMBOL), size = 1.0) +
  scale_x_continuous("False Positive Rate", limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous("True Positive Rate", limits = c(0,1), expand = c(0,0)) +
  ggtitle("ROC Curves for Selected Vitamin A Genes") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# Build legend labels with AUC
auc_labs <- auc_df |>
  mutate(label = sprintf("%s (AUC = %.2f)", SYMBOL, AUC)) |>
  select(SYMBOL, label)

p <- p +
  scale_color_discrete(labels = setNames(auc_labs$label, auc_labs$SYMBOL)) +
  scale_fill_discrete(labels = setNames(auc_labs$label, auc_labs$SYMBOL))

ggsave(out_png, p, width = 10, height = 6, dpi = 160)

# ---- Write AUC table ----
write_csv(auc_df, out_csv)
