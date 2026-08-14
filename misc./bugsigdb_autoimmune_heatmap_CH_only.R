# bugsigdb_autoimmune_LOR_heatmap.R
# ------------------------------------------------------------
# Microbe × Autoimmune condition signed log-odds ratio (LOR) heatmap
# - Robust ontology subsetting (EFO) for autoimmune disease (+ descendants)
# - GENUS level; maps specific taxa up to genus
# - Ensures ≥2 conditions, boosts variation, saves PNG + CSVs
#   * bugsigdb_autoimmune_LOR_matrix.csv        (raw LOR)
#   * bugsigdb_autoimmune_LOR_matrix_used.csv   (raw/centered LOR actually plotted)
#   * bugsigdb_autoimmune_LOR_heatmap.png       (figure with mid-tone emphasis)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(bugsigdbr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ontologyIndex)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(matrixStats)
  library(readr)
})

# -----------------------------
# Parameters (tune as needed)
# -----------------------------
PARAMS <- list(
  restrict_to_human       = TRUE,   # set FALSE to broaden
  restrict_to_feces       = FALSE,  # keep FALSE to broaden
  min_cond_signatures     = 2,      # min signatures per condition (auto-relaxes to ensure ≥2)
  min_taxon_occurrences   = 3,      # keep taxa seen in ≥3 autoimmune signatures overall
  tax_level               = "genus",
  pseudocount             = 0.5,
  top_taxa                = 100,    # set to Inf to keep all taxa
  use_devel               = FALSE,  # TRUE -> use "devel" export (bigger, unreviewed)
  # Plot choice: "raw" (recommended) or "row_centered"
  plot_mode               = "raw",
  # Quantile anchors for color range (symmetric around 0 using max(|q05|,|q95|))
  quantiles               = c(0.05, 0.95)
)

# -----------------------------
# Helpers
# -----------------------------
nz <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }

save_png <- function(file, width = 2400, height = 1800, res = 300) {
  sys <- tolower(Sys.info()[["sysname"]])
  if (sys == "windows") {
    png(file, width = width, height = height, res = res, type = "windows")
  } else {
    png(file, width = width, height = height, res = res)
  }
}

safe_stop <- function(msg) { stop(msg, call. = FALSE) }

# -----------------------------
# 1) Import BugSigDB (stable by default)
# -----------------------------
load_bug_sigdb <- function(use_devel = FALSE) {
  if (!use_devel) {
    df <- importBugSigDB()  # latest reviewed stable snapshot
  } else {
    df <- importBugSigDB(version = "devel", cache = FALSE)
  }
  if (nrow(df) == 0L) safe_stop("BugSigDB import returned 0 rows.")
  df
}

# -----------------------------
# 2) Autoimmune subset (Condition ~ EFO 'autoimmune disease' + descendants)
# -----------------------------
subset_autoimmune <- function(df) {
  efo <- getOntology("efo")
  idx_exact <- which(tolower(efo$name) == "autoimmune disease")
  if (length(idx_exact) > 0) {
    autoimmune_label <- unname(efo$name[idx_exact[1]])
  } else {
    ix <- grep("autoimmun", efo$name, ignore.case = TRUE)
    if (!length(ix)) safe_stop("Could not find an 'autoimmune' term in EFO.")
    autoimmune_label <- unname(efo$name[ix[1]])
  }
  
  df_ai <- tryCatch(
    subsetByOntology(df, column = "Condition", term = autoimmune_label, ontology = efo),
    error = function(e) NULL
  )
  if (is.null(df_ai) || nrow(df_ai) == 0L) {
    ai_id <- names(efo$name)[match(autoimmune_label, efo$name)]
    desc_ids <- unique(c(ai_id, ontologyIndex::get_descendants(efo, ai_id)))
    desc_names <- unname(efo$name[desc_ids])
    df_ai <- df[df$Condition %in% desc_names, , drop = FALSE]
  }
  if (nrow(df_ai) == 0L) safe_stop("Autoimmune ontology subset is empty.")
  df_ai
}

# -----------------------------
# 3) Filters (human; body-site) and ensure ≥2 conditions
# -----------------------------
filter_and_require_conditions <- function(df_ai, restrict_to_human, restrict_to_feces, min_cond_signatures) {
  if (restrict_to_human) {
    human_regex <- "(?i)^homo sapiens$|^human$"
    df_ai <- df_ai[grepl(human_regex, nz(df_ai$`Host species`)), , drop = FALSE]
  }
  if (restrict_to_feces) {
    feces_regex <- "(?i)feces|faeces|stool|fecal"
    df_ai <- df_ai[grepl(feces_regex, nz(df_ai$`Body site`)), , drop = FALSE]
  }
  if (nrow(df_ai) == 0L) safe_stop("All rows removed by human/body-site filters.")
  
  df_ai <- as.data.frame(df_ai, stringsAsFactors = FALSE)
  cond_counts <- df_ai %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(n_sig = dplyr::n(), .groups = "drop")
  
  keep_conds <- cond_counts$Condition[cond_counts$n_sig >= min_cond_signatures]
  if (length(keep_conds) < 2) {
    message(sprintf("[WARN] < 2 conditions meet min_cond_signatures = %d. Relaxing to 1.",
                    min_cond_signatures))
    keep_conds <- cond_counts$Condition[cond_counts$n_sig >= 1]
  }
  if (length(keep_conds) < 2) {
    message("[WARN] Still < 2 conditions after relaxing. Keeping the top two by signature count.")
    keep_conds <- utils::head(cond_counts$Condition[order(-cond_counts$n_sig)], 2)
  }
  
  df_ai <- df_ai[df_ai$Condition %in% keep_conds, , drop = FALSE]
  if (nrow(df_ai) == 0L || length(unique(df_ai$Condition)) < 2L) {
    safe_stop("After filtering, < 2 autoimmune conditions remain. Consider broadening filters or using 'devel'.")
  }
  df_ai
}

# -----------------------------
# 4) Build robust long table (Taxon, Condition, Dir) at GENUS level
# -----------------------------
build_long <- function(df_ai, tax_level) {
  sigs <- getSignatures(
    df_ai,
    tax.id.type     = "taxname",
    tax.level       = tax_level,
    exact.tax.level = FALSE,
    min.size        = 1
  )
  dir_sign <- ifelse(grepl("(?i)increase", nz(df_ai$`Abundance in Group 1`)),  1L,
                     ifelse(grepl("(?i)decrease", nz(df_ai$`Abundance in Group 1`)), -1L, 0L))
  lens  <- lengths(sigs)
  valid <- which(dir_sign != 0L & lens > 0L & !is.na(df_ai$Condition))
  
  df_ai_valid <- df_ai[valid, , drop = FALSE]
  dir_valid   <- dir_sign[valid]
  sigs_valid  <- sigs[valid]
  if (!length(sigs_valid)) safe_stop("No signatures left after direction/taxon filtering.")
  
  out <- vector("list", length(sigs_valid))
  k <- 0L
  for (i in seq_along(sigs_valid)) {
    tx <- sigs_valid[[i]]
    if (!length(tx) || is.na(dir_valid[i]) || dir_valid[i] == 0L || is.na(df_ai_valid$Condition[i])) next
    k <- k + 1L
    out[[k]] <- data.frame(
      Taxon     = as.character(tx),
      Condition = rep(df_ai_valid$Condition[i], length(tx)),
      Dir       = rep(dir_valid[i],            length(tx)),
      stringsAsFactors = FALSE
    )
  }
  if (k == 0L) safe_stop("After filtering for direction and mapped taxa, no rows remain.")
  long <- dplyr::bind_rows(out[seq_len(k)])
  list(long = long, df_ai = df_ai_valid)
}

# -----------------------------
# 5) Compute LOR matrix (Taxon × Condition)
# -----------------------------
lor_matrix <- function(long, min_taxon_occurrences, pseudocount) {
  counts <- long %>%
    dplyr::mutate(is_up = Dir > 0L, is_down = Dir < 0L) %>%
    dplyr::group_by(Taxon, Condition) %>%
    dplyr::summarise(UP = sum(is_up), DOWN = sum(is_down), .groups = "drop")
  
  totals_by_taxon <- counts %>%
    dplyr::group_by(Taxon) %>%
    dplyr::summarise(UP_all = sum(UP), DOWN_all = sum(DOWN), .groups = "drop")
  
  tot_occ   <- rowSums(as.matrix(totals_by_taxon[, c("UP_all","DOWN_all")]))
  keep_taxa <- totals_by_taxon$Taxon[tot_occ >= min_taxon_occurrences]
  counts          <- counts %>% dplyr::filter(Taxon %in% keep_taxa)
  totals_by_taxon <- totals_by_taxon %>% dplyr::filter(Taxon %in% keep_taxa)
  if (!nrow(counts)) safe_stop("No taxa left after 'min_taxon_occurrences' filtering.")
  
  counts <- counts %>%
    dplyr::left_join(totals_by_taxon, by = "Taxon") %>%
    dplyr::mutate(
      UP_others   = pmax(UP_all   - UP,   0L),
      DOWN_others = pmax(DOWN_all - DOWN, 0L),
      LOR = log( ((UP + pseudocount) * (DOWN_others + pseudocount)) /
                   ((DOWN + pseudocount) * (UP_others + pseudocount)) )
    )
  
  LOR_wide <- counts %>%
    dplyr::select(Taxon, Condition, LOR) %>%
    tidyr::pivot_wider(names_from = Condition, values_from = LOR, values_fill = 0) %>%
    as.data.frame()
  rownames(LOR_wide) <- LOR_wide$Taxon
  LOR_wide$Taxon <- NULL
  M <- as.matrix(LOR_wide)
  
  if (ncol(M) > 1) M <- M[, colSums(abs(M), na.rm = TRUE) > 0, drop = FALSE]
  if (nrow(M) == 0L || ncol(M) == 0L) safe_stop("LOR matrix is empty after cleaning.")
  M
}

# -----------------------------
# 6) Plot heatmap with explicit 0 color (mid-tone emphasis, no row z-scaling)
# -----------------------------
plot_LOR_midtones <- function(
    M, df_ai, outfile_png, top_taxa,
    plot_mode = c("raw","row_centered"),
    probs = c(0.05, 0.95),
    zero_col = "#F2F2F2",            # <-- explicit neutral color for 0
    cell_border = "#DADADA"          # subtle gridlines so zero cells are visible
) {
  plot_mode <- match.arg(plot_mode)
  
  # Option 1: raw LOR (recommended). Option 2: row-centered (subtract row means).
  if (plot_mode == "raw") {
    M_plot <- M
  } else {
    M_plot <- M - rowMeans(M, na.rm = TRUE)
  }
  
  # Keep top taxa by total |raw LOR| for readability (rank by raw magnitude)
  if (is.finite(top_taxa) && nrow(M_plot) > top_taxa) {
    total_abs <- rowSums(abs(M), na.rm = TRUE)
    ord <- order(total_abs, decreasing = TRUE)
    M_plot <- M_plot[ord[seq_len(top_taxa)], , drop = FALSE]
  }
  
  # Data-driven symmetric color limits using quantiles
  vals <- as.numeric(M_plot); vals <- vals[is.finite(vals)]
  if (!length(vals)) stop("No finite LOR values to plot.")
  q <- stats::quantile(vals, probs = probs, na.rm = TRUE)
  L <- max(abs(q))
  if (!is.finite(L) || L == 0) {
    L <- max(abs(vals), na.rm = TRUE)
    if (!is.finite(L) || L == 0) L <- 1
  }
  
  # Mid-tone emphasized 5-stop palette with explicit neutral at 0
  col_fun <- circlize::colorRamp2(
    c(-L, -L*0.33, 0, L*0.33, L),
    c("#053061", "#2166AC", zero_col, "#B2182B", "#67001F")
  )
  
  # Study-count bar
  study_counts <- df_ai %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(n_studies = dplyr::n_distinct(.data$Study), .groups = "drop")
  sc_vec <- setNames(study_counts$n_studies, study_counts$Condition)
  study_vec <- sc_vec[colnames(M_plot)]
  study_vec[is.na(study_vec)] <- 0
  
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Studies = ComplexHeatmap::anno_barplot(study_vec, gp = grid::gpar(fill = "#6C757D", col = NA),
                                           border = FALSE, height = grid::unit(2.0, "cm")),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontface = "bold")
  )
  
  ht <- ComplexHeatmap::Heatmap(
    M_plot, name = "LOR", col = col_fun,
    na_col = "#FFFFFF",                       # NAs are white; 0 is light gray (zero_col)
    rect_gp = grid::gpar(col = cell_border, lwd = 0.3),  # subtle gridlines
    cluster_rows = TRUE, cluster_columns = TRUE,
    row_title = "Microbes (genus)",
    column_title = "Autoimmune conditions (bar = # of studies)",
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = grid::gpar(cex = 0.6), column_names_gp = grid::gpar(cex = 0.7),
    top_annotation = ha,
    heatmap_legend_param = list(
      title = "LOR",
      at = c(-L, -L*0.5, 0, L*0.5, L),         # explicit 0 tick
      labels = sprintf(c("-%0.2f","-%0.2f","0","+%0.2f","+%0.2f"), L, L*0.5, L*0.5, L),
      legend_direction = "horizontal", title_position = "topcenter",
      grid_width = grid::unit(5, "mm"), grid_height = grid::unit(5, "mm")
    )
  )
  
  # Save PNG (portable—no Cairo requirement)
  save_png(outfile_png)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  dev.off()
  
  # Save matrix actually plotted (for provenance)
  readr::write_csv(
    as.data.frame(M_plot) %>% tibble::rownames_to_column("Taxon"),
    "bugsigdb_autoimmune_LOR_matrix_used.csv"
  )
  
  message(sprintf("✓ LOR heatmap saved: %s", normalizePath(outfile_png, winslash = "/")))
  message(sprintf("[INFO] Using %d autoimmune conditions: %s",
                  ncol(M_plot), paste(colnames(M_plot), collapse = " | ")))
  invisible(TRUE)
}

# Run if executed (not just sourced)
if (sys.nframe() == 0) {
  main()
}