# ==========================================
# Autoimmune – Non‑Immune Biomarkers (End-to-End; curatedMetagenomicData 3.14.0 compatible)
# ==========================================

options(warn = 1)  # warnings should not abort

suppressPackageStartupMessages({
  library(rentrez)   # NCBI E-utilities in R (programmatic BioProject/SRA)  [docs]  [1](https://andersvercelli.com/packages/release/data/experiment/manuals/curatedMetagenomicData/man/curatedMetagenomicData.pdf)
  library(httr)
  library(readr)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(glmnet)    # elastic-net  [6](https://bioc.r-universe.dev/phyloseq/doc/manual.html)
  library(pROC)      # ROC/AUC      [7](https://cran.r-project.org/web//packages/glmnet/glmnet.pdf)
})

# ---------- checks for curatedMetagenomicData ----------
have_cmd <- requireNamespace("curatedMetagenomicData", quietly = TRUE) &&
  requireNamespace("SummarizedExperiment", quietly = TRUE)
if (!have_cmd) {
  message("NOTE: curatedMetagenomicData/SummarizedExperiment not found. ",
          "RunInfo CSVs will be produced; features/modeling skipped. ",
          "Install via BiocManager::install('curatedMetagenomicData').")
}

# ---------- small utilities ----------
safe_write_csv <- function(df, path) {
  if (!"Symbol" %in% names(df)) df <- df %>% mutate(Symbol = .data[[1]]) %>% relocate(Symbol)
  readr::write_csv(df, path)
  message("Wrote: ", normalizePath(path, winslash = "/"))
}

bind_runinfo <- function(lst) {
  lst <- purrr::keep(lst, ~ is.data.frame(.x) && nrow(.x) > 0)
  if (!length(lst)) return(tibble(Run = character()))
  lst <- purrr::map(lst, ~ dplyr::mutate(.x, dplyr::across(dplyr::everything(), as.character)))
  purrr::list_rbind(lst)
}

standardize_run_column <- function(df, bp_accession = NA_character_) {
  if (!"Run" %in% names(df)) {
    cand <- which(grepl("^run(_)?accession$|^run$", names(df), ignore.case = TRUE))
    if (length(cand) == 1) names(df)[cand] <- "Run"
  }
  if (!"Run" %in% names(df)) {
    warning("No recognizable Run column for BioProject: ", bp_accession)
    return(NULL)
  }
  df
}

dedup_by <- function(df, key) {
  if (key %in% names(df)) df[!duplicated(df[[key]]), , drop = FALSE] else df
}

# ---------- SRA RunInfo fetch (Run Selector + E-utilities fallback) ----------
fetch_runinfo_for_bioproject <- function(bp_accession, pause_sec = 0.34) {
  # Mirrors the SRA Run Selector "RunInfo Table" download; falls back to E-utilities  [3](https://biopython.org/docs/dev/Tutorial/chapter_entrez.html)
  base <- "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi"
  url  <- paste0(base, "?save=efetch&db=sra&rettype=runinfo&term=", bp_accession)
  
  get_and_parse <- function(u) {
    resp <- try(httr::GET(u), silent = TRUE)
    if (inherits(resp, "try-error") || httr::http_error(resp)) return(NULL)
    txt  <- httr::content(resp, as = "text", encoding = "UTF-8")
    if (is.null(txt) || !nzchar(txt)) return(NULL)
    suppressMessages(readr::read_csv(txt, col_types = readr::cols(.default = readr::col_character())))
  }
  
  df <- get_and_parse(url)
  
  if (is.null(df) || !nrow(df)) {
    # E-utilities fallback via rentrez  [2](https://deepwiki.com/waldronlab/curatedMetagenomicData/2.3-sample-selection-with-returnsamples)
    s_ids <- try(rentrez::entrez_search(db = "sra",
                                        term = paste0(bp_accession, "[bioproject]"),
                                        retmax = 999999),
                 silent = TRUE)
    if (!inherits(s_ids, "try-error") && length(s_ids$ids)) {
      csv_txt <- try(rentrez::entrez_fetch(db = "sra",
                                           id = s_ids$ids,
                                           rettype = "runinfo",
                                           retmode = "text"),
                     silent = TRUE)
      if (!inherits(csv_txt, "try-error") && is.character(csv_txt) && nzchar(csv_txt)) {
        df <- suppressMessages(readr::read_csv(csv_txt, col_types = readr::cols(.default = readr::col_character())))
      }
    }
  }
  
  if (is.null(df) || !nrow(df)) return(NULL)
  df <- standardize_run_column(df, bp_accession)
  if (is.null(df)) return(NULL)
  df$BioProject <- bp_accession
  Sys.sleep(pause_sec)
  df
}

# ---------- helpers for curatedMetagenomicData ----------
find_bioproject_col <- function(df) {
  # Robust PRJNA/PRJEB/PRJDB detection; returns NA if none
  prj_regex <- "PRJ[A-Z]{1,4}\\d+"
  nm <- names(df)
  is_chr <- vapply(nm, function(n) is.character(df[[n]]) || is.factor(df[[n]]), logical(1))
  nm_chr <- nm[is_chr]
  if (!length(nm_chr)) return(NA_character_)
  has_prj <- vapply(
    nm_chr,
    function(n) any(grepl(prj_regex, as.character(df[[n]]), ignore.case = TRUE), na.rm = TRUE),
    logical(1)
  )
  if (any(has_prj)) nm_chr[which(has_prj)[1]] else NA_character_
}

count_dplyr <- function(.data, ...) dplyr::count(.data, ...)

find_col <- function(df, patterns) {
  nm <- names(df)
  i  <- which(vapply(patterns, function(p) any(grepl(p, nm, ignore.case = TRUE)), logical(1)))
  if (length(i) == 0) return(NA_character_)
  pat <- patterns[i[1]]
  hit <- nm[grepl(pat, nm, ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
}

# ==========================================
# 1) SRA — reproduce two BioProject queries + RunInfo (E-utilities)
# ==========================================
term_microbiome <- '("autoimmune"[All Fields] AND microbiome[All Fields]) AND bioproject_sra[filter] NOT bioproject_gap[filter]'
term_diet       <- '("autoimmune"[All Fields] AND diet[All Fields]) AND bioproject_sra[filter] NOT bioproject_gap[filter]'

bp_search <- function(term) {
  s <- entrez_search(db = "bioproject", term = term, retmax = 99999, use_history = TRUE)
  tibble(bioproject_uid = s$ids, count = s$count,
         webenv = s$web_history$WebEnv, query_key = s$web_history$QueryKey)
}

bp_micro <- bp_search(term_microbiome)   # rentrez: E-utilities in R  [1](https://andersvercelli.com/packages/release/data/experiment/manuals/curatedMetagenomicData/man/curatedMetagenomicData.pdf)
bp_diet  <- bp_search(term_diet)

bp_summary_to_accession <- function(uids) {
  if (!length(uids)) return(character(0))
  sums <- entrez_summary(db = "bioproject", id = uids)
  accs <- purrr::map_chr(sums, ~ .x$project_acc)
  unique(accs)
}

bp_micro_acc <- bp_summary_to_accession(bp_micro$bioproject_uid)
bp_diet_acc  <- bp_summary_to_accession(bp_diet$bioproject_uid)

runinfo_micro_list <- purrr::map(bp_micro_acc, fetch_runinfo_for_bioproject)
runinfo_diet_list  <- purrr::map(bp_diet_acc,  fetch_runinfo_for_bioproject)

runinfo_micro <- bind_runinfo(runinfo_micro_list)
runinfo_diet  <- bind_runinfo(runinfo_diet_list)

# Base‑R de‑dup (no any_of()), sorted by Run
if (nrow(runinfo_micro) == 0) {
  warning("No RunInfo rows were retrieved for 'autoimmune AND microbiome'.")
} else {
  runinfo_micro <- dedup_by(runinfo_micro, "Run")
  runinfo_micro <- runinfo_micro[order(runinfo_micro[["Run"]]), , drop = FALSE]
  safe_write_csv(dplyr::mutate(runinfo_micro, Symbol = Run), "sra_runinfo_autoimmune_microbiome.csv")
}

if (nrow(runinfo_diet) == 0) {
  warning("No RunInfo rows were retrieved for 'autoimmune AND diet'.")
} else {
  runinfo_diet <- dedup_by(runinfo_diet, "Run")
  runinfo_diet <- runinfo_diet[order(runinfo_diet[["Run"]]), , drop = FALSE]
  safe_write_csv(dplyr::mutate(runinfo_diet, Symbol = Run), "sra_runinfo_autoimmune_diet.csv")
}

# ==========================================
# 2) curatedMetagenomicData — sample-based retrieval via returnSamples()
# ==========================================
if (have_cmd) {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  
  smd_sym <- try(getExportedValue("curatedMetagenomicData", "sampleMetadata"), silent = TRUE)
  if (inherits(smd_sym, "try-error")) stop("curatedMetagenomicData::sampleMetadata not found.")
  smd <- if (is.function(smd_sym)) smd_sym() else smd_sym
  smd <- tibble::as_tibble(smd)
  
  autoimmune_terms <- paste(c("autoimmune","IBD","inflammatory bowel","Crohn",
                              "ulcerative colitis","rheumatoid","multiple sclerosis",
                              "systemic lupus","psoriasis","type 1 diabetes","celiac"),
                            collapse = "|")
  diet_terms <- paste(c("diet","dietary","high[- ]fiber","FODMAP","keto","vegan",
                        "vegetarian","Mediterranean","western","polyphenol","prebiotic","probiotic","food"),
                      collapse = "|")
  
  smd_chr  <- smd %>% mutate(across(.cols = everything(), ~ as.character(.)))
  smd_blob <- smd_chr %>% unite(".search_blob", everything(), sep = " ", remove = FALSE)
  bp_col   <- find_bioproject_col(smd_blob)
  
  smd_blob <- smd_blob %>%
    mutate(
      match_autoimmune = str_detect(tolower(.search_blob), tolower(autoimmune_terms)),
      match_diet       = str_detect(tolower(.search_blob), tolower(diet_terms)),
      in_our_sra       = if (!is.na(bp_col)) .data[[bp_col]] %in% c(bp_micro_acc, bp_diet_acc) else FALSE
    )
  
  candidate_cmd <- smd_blob %>%
    filter(match_autoimmune | match_diet | in_our_sra) %>%
    tibble::as_tibble()
  
  print(count_dplyr(candidate_cmd, match_autoimmune, match_diet, in_our_sra))
  
  # returnSamples() requires 'study_name' + 'sample_id' (no 'study=' arg)  [4](https://www.metafor-project.org/doku.php/help)
  study_col <- if ("study_name" %in% names(candidate_cmd)) {
    "study_name"
  } else {
    find_col(candidate_cmd, c("^study_name$", "^study$", "dataset"))
  }
  
  sample_col <- if ("sample_id" %in% names(candidate_cmd)) {
    "sample_id"
  } else {
    find_col(candidate_cmd, c("^sample_id$", "^sample$", "^id$"))
  }
  
  if (is.na(study_col) || is.na(sample_col)) {
    warning("Could not find 'study_name' and/or 'sample_id' in candidate metadata; skipping cMD retrieval.")
    have_cmd <- FALSE
  }
  
  if (have_cmd) {
    candidate_cmd_clean <- candidate_cmd %>%
      select(any_of(c(study_col, sample_col)), where(~ !all(is.na(.x)))) %>%
      filter(!is.na(.data[[study_col]]), !is.na(.data[[sample_col]])) %>%
      rename(study_name = !!study_col, sample_id = !!sample_col) %>%
      distinct(study_name, sample_id, .keep_all = TRUE) %>%
      as_tibble()
    
    # NA cleanup without across(); do not touch keys
    num_cols <- setdiff(names(candidate_cmd_clean)[sapply(candidate_cmd_clean, is.numeric)],
                        c("study_name","sample_id"))
    chr_cols <- setdiff(names(candidate_cmd_clean)[sapply(candidate_cmd_clean, is.character)],
                        c("study_name","sample_id"))
    if (length(num_cols)) candidate_cmd_clean[num_cols] <- lapply(candidate_cmd_clean[num_cols], tidyr::replace_na, 0)
    if (length(chr_cols)) candidate_cmd_clean[chr_cols] <- lapply(candidate_cmd_clean[chr_cols], tidyr::replace_na, "")
    
    # ---- Temporarily silence package-internal lifecycle warning from returnSamples() ----
    old_lifecycle <- getOption("lifecycle_verbosity")
    options(lifecycle_verbosity = "quiet")
    
    # Sample-based retrieval by dataType (species & pathways)  [4](https://www.metafor-project.org/doku.php/help)
    species_list <- curatedMetagenomicData::returnSamples(candidate_cmd_clean, dataType = "relative_abundance")
    pathway_list <- curatedMetagenomicData::returnSamples(candidate_cmd_clean, dataType = "pathway_abundance")
    
    # restore lifecycle verbosity
    options(lifecycle_verbosity = old_lifecycle)
    
    # ---------- normalize outputs to base lists and convert safely ----------
    to_base_list <- function(x) {
      if (is.null(x)) return(list())
      if (inherits(x, c("SimpleList","List"))) return(as.list(x))
      if (is.list(x)) return(x)
      list(x)
    }
    species_items <- to_base_list(species_list)
    pathway_items <- to_base_list(pathway_list)
    
    is_valid_se <- function(x) {
      tryCatch({
        inherits(x, c("SummarizedExperiment","RangedSummarizedExperiment","TreeSummarizedExperiment")) &&
          length(SummarizedExperiment::assayNames(x)) > 0 &&
          {
            an <- SummarizedExperiment::assayNames(x)[1]
            !is.na(an) &&
              !is.null(SummarizedExperiment::assay(x, an)) &&
              nrow(SummarizedExperiment::assay(x, an)) > 0 &&
              ncol(SummarizedExperiment::assay(x, an)) > 0
          }
      }, error = function(e) FALSE)
    }
    
    to_long_safe <- function(se_obj, feature_type) {
      tryCatch({
        an <- SummarizedExperiment::assayNames(se_obj)[1]
        A  <- SummarizedExperiment::assay(se_obj, an)
        meta <- SummarizedExperiment::colData(se_obj) |>
          as.data.frame() |>
          tibble::rownames_to_column("Sample")
        tibble::as_tibble(A, rownames = "Feature") |>
          tidyr::pivot_longer(-Feature, names_to = "Sample", values_to = "abundance") |>
          dplyr::left_join(meta, by = "Sample") |>
          dplyr::mutate(feature_type = feature_type)
      }, error = function(e) NULL)
    }
    
    as_long_tbl <- function(items, feature_type) {
      out <- vector("list", length(items)); k <- 0L
      for (i in seq_along(items)) {
        if (!is_valid_se(items[[i]])) next
        tmp <- to_long_safe(items[[i]], feature_type)
        if (!is.null(tmp) && nrow(tmp) > 0) { k <- k + 1L; out[[k]] <- tmp }
      }
      if (k == 0L) tibble() else dplyr::bind_rows(out[seq_len(k)])
    }
    
    species_tbls <- as_long_tbl(species_items, "species")
    pathway_tbls <- as_long_tbl(pathway_items, "pathway")
    
    message("#species rows: ", nrow(species_tbls), " | #pathway rows: ", nrow(pathway_tbls))
    if (nrow(species_tbls) == 0 && nrow(pathway_tbls) == 0) {
      warning("No species or pathway features returned—relax filters or inspect labels.")
      have_cmd <- FALSE
    }
    
    # ==========================================
    # 3) Combine, non-immune filter, matrices, LOSO modeling
    # ==========================================
    if (have_cmd) {
      clean_tbl <- dplyr::bind_rows(species_tbls, pathway_tbls)
      
      # optional: keep gut/stool
      if ("body_site" %in% names(clean_tbl)) {
        clean_tbl <- clean_tbl %>% filter(is.na(body_site) | grepl("stool|fec|gut", tolower(body_site)))
      }
      
      # infer Case (1 = autoimmune)
      guess_case <- function(df) {
        dx_cols <- c("disease","disease_subtype","phenotype","condition","group","diagnosis","label","status")
        have    <- dx_cols[dx_cols %in% names(df)]
        if (!length(have)) return(rep(NA_integer_, nrow(df)))
        dx <- df %>% unite(".dx", all_of(have), sep = " ", remove = FALSE)
        as.integer(str_detect(tolower(dx$.dx),
                              "autoimmune|ibd|crohn|ulcerative colitis|rheumatoid|multiple sclerosis|systemic lupus|psoriasis|type 1 diabetes|celiac"))
      }
      clean_tbl <- clean_tbl %>% mutate(Case = guess_case(.))
      
      safe_write_csv(clean_tbl %>% arrange(feature_type, Feature) %>% mutate(Symbol = Feature) %>% relocate(Symbol),
                     "cmd_features_long.csv")
      
      immune_regex <- "(immune|interferon|interleukin|antigen|chemokine|mhc|t[- ]cell|b[- ]cell|nk cell|complement)"
      species_wide <- clean_tbl %>% filter(feature_type == "species") %>%
        select(Symbol = Feature, Sample, abundance, Case)
      pathway_wide <- clean_tbl %>% filter(feature_type == "pathway") %>%
        filter(!str_detect(tolower(Feature), immune_regex)) %>%
        select(Symbol = Feature, Sample, abundance, Case)
      
      to_wide <- function(df) {
        if (!nrow(df)) return(list(mat = NULL, meta = NULL))
        mat <- df %>% select(Symbol, Sample, abundance) %>%
          pivot_wider(names_from = Symbol, values_from = abundance, values_fill = 0)
        meta <- df %>% distinct(Sample, Case)
        list(mat = mat, meta = meta)
      }
      
      sp <- to_wide(species_wide)
      pw <- to_wide(pathway_wide)
      
      if (is.null(sp$mat) || is.null(pw$mat)) {
        warning("Insufficient features to build matrices; modeling skipped.")
        have_cmd <- FALSE
      } else {
        common_samples <- intersect(sp$meta$Sample, pw$meta$Sample)
        X_species  <- sp$mat %>% filter(Sample %in% common_samples) %>% arrange(Sample)
        X_pathway  <- pw$mat %>% filter(Sample %in% common_samples) %>% arrange(Sample)
        meta_df    <- sp$meta %>% filter(Sample %in% common_samples) %>% arrange(Sample)
        
        stopifnot(identical(X_species$Sample, X_pathway$Sample),
                  identical(X_species$Sample, meta_df$Sample))
        
        X <- X_species %>% select(-Sample) %>% as.matrix()
        Y <- X_pathway %>% select(-Sample) %>% as.matrix()
        
        keep_cols <- function(M, thr = 0.10) {
          if (is.null(M) || nrow(M) == 0) return(M)
          keep <- colMeans(M > 0, na.rm = TRUE) >= thr
          if (!any(keep)) return(M[, 0, drop = FALSE])
          M[, keep, drop = FALSE]
        }
        X <- keep_cols(X, 0.10)
        Y <- keep_cols(Y, 0.10)
        
        if (ncol(X) + ncol(Y) < 2) {
          warning("Too few features after prevalence filtering; modeling skipped.")
          have_cmd <- FALSE
        } else {
          XY <- cbind(X, Y)
          XY_scaled <- scale(XY)
          
          study_col_guess <- c("study_name","study","dataset_name","Study")[c("study_name","study","dataset_name","Study") %in% names(clean_tbl)][1]
          if (!is.null(study_col_guess)) {
            study_lookup <- clean_tbl %>% distinct(Sample, Study = .data[[study_col_guess]])
            meta_df <- meta_df %>% left_join(study_lookup, by = "Sample")
          } else {
            meta_df$Study <- "ALL"
          }
          
          safe_write_csv(
            as_tibble(XY_scaled) %>% mutate(Sample = meta_df$Sample) %>% relocate(Sample) %>%
              mutate(Symbol = Sample) %>% relocate(Symbol) %>% arrange(Symbol),
            "combined_feature_matrix_scaled.csv"
          )
          safe_write_csv(meta_df %>% mutate(Symbol = Sample) %>% relocate(Symbol) %>% arrange(Symbol),
                         "combined_feature_metadata.csv")
          
          # ---------- LOSO elastic-net + AUC ----------
          set.seed(1)
          alpha_elastic <- 0.5
          y_bin <- meta_df$Case
          y_bin[is.na(y_bin)] <- 0L
          
          study_vec <- ifelse(is.na(meta_df$Study), "ALL", as.character(meta_df$Study))
          unique_studies <- unique(study_vec)
          pred_all <- rep(NA_real_, length(y_bin))
          coef_list <- list()
          
          for (left_out in unique_studies) {
            idx_test  <- which(study_vec == left_out)
            idx_train <- setdiff(seq_along(y_bin), idx_test)
            if (length(unique(y_bin[idx_train])) < 2) next
            
            cvfit <- cv.glmnet(x = XY_scaled[idx_train, , drop = FALSE], y = y_bin[idx_train],
                               family = "binomial", alpha = alpha_elastic, nfolds = 5, type.measure = "auc")
            fit   <- glmnet(x = XY_scaled[idx_train, , drop = FALSE], y = y_bin[idx_train],
                            family = "binomial", alpha = alpha_elastic, lambda = cvfit$lambda.min)
            
            pred_all[idx_test] <- as.numeric(predict(fit, newx = XY_scaled[idx_test, , drop = FALSE], type = "response"))
            
            co <- coef(fit); nz <- which(co != 0)
            coef_df <- tibble(
              Feature = rownames(co)[nz],
              Coef    = as.numeric(co[nz]),
              Study_left_out = left_out
            ) %>% filter(Feature != "(Intercept)")
            coef_list[[left_out]] <- coef_df
          }
          
          roc_obj <- pROC::roc(response = y_bin, predictor = pred_all, quiet = TRUE)
          auc_val <- as.numeric(pROC::auc(roc_obj))
          message(sprintf("Overall LOSO AUC = %.3f", auc_val))
          
          study_perf <- tibble(Study = unique_studies) %>%
            rowwise() %>%
            mutate(AUC = {
              idx <- which(study_vec == Study)
              if (length(unique(y_bin[idx])) < 2 || all(is.na(pred_all[idx]))) NA_real_
              else as.numeric(pROC::auc(pROC::roc(y_bin[idx], pred_all[idx], quiet = TRUE)))
            }) %>% ungroup()
          
          safe_write_csv(study_perf %>% mutate(Symbol = Study) %>% relocate(Symbol) %>% arrange(Symbol),
                         "validation_LOSO_by_study_auc.csv")
          
          coef_tbl <- bind_rows(coef_list) %>%
            group_by(Feature) %>%
            summarize(n_selected = n(), mean_coef = mean(Coef), .groups = "drop") %>%
            arrange(desc(n_selected), desc(abs(mean_coef)))
          
          safe_write_csv(coef_tbl %>% mutate(Symbol = Feature) %>% relocate(Symbol) %>% arrange(Symbol),
                         "biomarkers_nonimmune_elasticnet.csv")
          
          cv_all <- cv.glmnet(x = XY_scaled, y = y_bin, family = "binomial",
                              alpha = alpha_elastic, nfolds = 5, type.measure = "auc")
          final_fit <- glmnet(x = XY_scaled, y = y_bin, family = "binomial",
                              alpha = alpha_elastic, lambda = cv_all$lambda.min)
          
          saveRDS(list(model = final_fit,
                       columns = colnames(XY_scaled),
                       center = attr(XY_scaled, "scaled:center"),
                       scale  = attr(XY_scaled, "scaled:scale")),
                  file = "Autoimmune_NonImmune_BiomarkerScore_model.rds")
          message("Saved model -> Autoimmune_NonImmune_BiomarkerScore_model.rds")
        }
      }
    }
  }
}

message("Pipeline finished.")
``

# ==========================================
# Compare SRA metadata classifications between two RunInfo lists
# Outputs: sra_metadata_classifications_comparison.csv
# ==========================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

# -------------------------------
# 0) Helper: fail early if files missing
# -------------------------------
expect_file <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
}

# -------------------------------
# 1) Read both RunInfo CSVs (force character types to avoid clashes)
# -------------------------------
read_runinfo <- function(path, cohort_label) {
  expect_file(path)
  df <- suppressMessages(readr::read_csv(path, col_types = readr::cols(.default = readr::col_character())))
  # If a previous step already added Symbol, keep it; doesn't hurt
  df$cohort <- cohort_label
  df
}

micro <- read_runinfo("sra_runinfo_autoimmune_microbiome.csv", "autoimmune_microbiome")
diet  <- read_runinfo("sra_runinfo_autoimmune_diet.csv",       "autoimmune_diet")

# Combine for unified handling
ri <- dplyr::bind_rows(micro, diet)

# -------------------------------
# 2) Normalize column names → canonical SRA labels
#    (Handles casing and common alternates)
# -------------------------------
normalize_colnames <- function(nm) {
  map <- c(
    "librarystrategy"   = "LibraryStrategy",
    "librarysource"     = "LibrarySource",
    "libraryselection"  = "LibrarySelection",
    "librarylayout"     = "LibraryLayout",
    "platform"          = "Platform",
    "model"             = "Model",
    "centername"        = "CenterName",
    "scientificname"    = "ScientificName",
    "srastudy"          = "SRAStudy",
    "bioproject"        = "BioProject",
    "biosample"         = "BioSample",
    "samplename"        = "SampleName",
    "run"               = "Run",
    "experiment"        = "Experiment",
    "instrument"        = "Model"          # sometimes appears as "Instrument"
  )
  std <- nm
  lower <- tolower(nm)
  for (i in seq_along(nm)) {
    if (lower[i] %in% names(map)) std[i] <- map[[lower[i]]]
  }
  std
}
names(ri) <- normalize_colnames(names(ri))

# -------------------------------
# 3) Choose metadata classification variables to compare
#    (will keep only those that exist in your files)
# -------------------------------
candidate_class_vars <- c(
  "LibraryStrategy", "LibrarySource", "LibrarySelection", "LibraryLayout",
  "Platform", "Model", "CenterName", "ScientificName", "SRAStudy", "BioProject"
)

class_vars <- intersect(candidate_class_vars, names(ri))
if (!length(class_vars)) {
  stop("No expected classification columns were found. Available columns are: ",
       paste(names(ri), collapse = ", "), call. = FALSE)
}

# Keep only cohort + the classification variables
ri_sel <- ri %>% dplyr::select(cohort, dplyr::all_of(class_vars)) %>% tibble::as_tibble()

# -------------------------------
# 4) Clean obvious NA-like strings to real NA
# -------------------------------
clean_na <- function(x) {
  x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
  x
}
ri_sel[class_vars] <- lapply(ri_sel[class_vars], clean_na)

# -------------------------------
# 5) Long-format counts per cohort (masking-safe dplyr::count)
# -------------------------------
long_df <- ri_sel %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(class_vars),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  dplyr::filter(!is.na(value))

# Use dplyr::count explicitly to avoid plyr masking
long_counts <- long_df %>%
  dplyr::count(variable, value, cohort, name = "n", sort = FALSE) %>%
  tidyr::pivot_wider(
    names_from  = cohort,
    values_from = n,
    values_fill = 0
  )

# Ensure both cohorts are present as columns (even if one is absent)
for (cn in c("autoimmune_microbiome","autoimmune_diet")) {
  if (!(cn %in% names(long_counts))) long_counts[[cn]] <- 0L
}

# -------------------------------
# 6) Add totals & percentages within each variable
# -------------------------------
totals <- long_counts %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    total_microbiome = sum(.data$autoimmune_microbiome, na.rm = TRUE),
    total_diet       = sum(.data$autoimmune_diet,       na.rm = TRUE),
    .groups          = "drop"
  )

summary_tbl <- long_counts %>%
  dplyr::left_join(totals, by = "variable") %>%
  dplyr::mutate(
    microbiome_pct = ifelse(total_microbiome > 0,
                            100 * autoimmune_microbiome / total_microbiome, 0),
    diet_pct       = ifelse(total_diet > 0,
                            100 * autoimmune_diet       / total_diet,       0),
    total          = autoimmune_microbiome + autoimmune_diet
  ) %>%
  # Put Symbol first to match your preferred CSV layout
  dplyr::transmute(
    Symbol          = value,
    variable,
    autoimmune_microbiome,
    microbiome_pct  = round(microbiome_pct, 2),
    autoimmune_diet,
    diet_pct        = round(diet_pct, 2),
    total
  ) %>%
  dplyr::arrange(variable, dplyr::desc(total), Symbol)

# -------------------------------
# 7) Write the comparison table + console preview
# -------------------------------
out_path <- "sra_metadata_classifications_comparison.csv"
readr::write_csv(summary_tbl, out_path)
message("Wrote: ", normalizePath(out_path, winslash = "/"))

# Preview: top 10 values per variable by total
preview <- summary_tbl %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(rank = dplyr::row_number(dplyr::desc(total))) %>%
  dplyr::filter(rank <= 10) %>%
  dplyr::ungroup() %>%
  dplyr::select(-rank)

print(preview, n = 60, width = 120)
``

# -------------------------------
# 1) Read SRA RunInfo (all-character) and combine
# -------------------------------
read_runinfo <- function(path, cohort_label) {
  stopifnot(file.exists(path))
  df <- suppressMessages(readr::read_csv(path, col_types = readr::cols(.default = readr::col_character())))
  df$cohort <- cohort_label
  df
}

micro <- read_runinfo("sra_runinfo_autoimmune_microbiome.csv", "autoimmune_microbiome")
diet  <- read_runinfo("sra_runinfo_autoimmune_diet.csv",       "autoimmune_diet")
ri    <- dplyr::bind_rows(micro, diet)

# -------------------------------
# 2) Normalize column names to canonical SRA labels
# -------------------------------
normalize_colnames <- function(nm) {
  map <- c(
    "librarystrategy"="LibraryStrategy","librarysource"="LibrarySource","libraryselection"="LibrarySelection",
    "librarylayout"="LibraryLayout","platform"="Platform","model"="Model","instrument"="Model",
    "centername"="CenterName","scientificname"="ScientificName","srastudy"="SRAStudy","bioproject"="BioProject",
    "biosample"="BioSample","samplename"="SampleName","run"="Run","experiment"="Experiment",
    "avglength"="avgLength"
  )
  std <- nm; lower <- tolower(nm)
  for (i in seq_along(nm)) if (lower[i] %in% names(map)) std[i] <- map[[lower[i]]]
  std
}
names(ri) <- normalize_colnames(names(ri))

# -------------------------------
# 3) Build the "most important" classification set
#    + derived ReadLengthBin if avgLength exists
# -------------------------------
important_vars <- c("LibraryStrategy","LibraryLayout","Platform","Model")
has_avglen <- "avgLength" %in% names(ri)

# Safe numeric parse for avgLength; build bins if present
if (has_avglen) {
  # parse_number handles "150", "150.0", "150 bp", etc.
  rl <- suppressWarnings(readr::parse_number(ri$avgLength))
  # Construct bins that are informative for metagenomics profiling
  rl_bin <- cut(
    rl,
    breaks = c(-Inf, 75, 150, 250, 500, Inf),
    labels = c("≤75", "76–150", "151–250", "251–500", ">500"),
    right = TRUE
  )
  ri$ReadLengthBin <- as.character(rl_bin)
  important_vars <- c(important_vars, "ReadLengthBin")
}

# Keep only cohort + important fields that actually exist
keep_vars <- intersect(important_vars, names(ri))
stopifnot(length(keep_vars) > 0)
ri_sel <- ri %>% dplyr::select(cohort, dplyr::all_of(keep_vars)) %>% tibble::as_tibble()

# Clean obvious NA-like strings
clean_na <- function(x) { x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_; x }
ri_sel[keep_vars] <- lapply(ri_sel[keep_vars], clean_na)

# -------------------------------
# 4) Long counts (masking-safe) and cohort totals
# -------------------------------
long_df <- ri_sel %>%
  tidyr::pivot_longer(cols = dplyr::all_of(keep_vars),
                      names_to = "variable", values_to = "value") %>%
  dplyr::filter(!is.na(value))

long_counts <- long_df %>%
  dplyr::count(variable, value, cohort, name = "n") %>%
  tidyr::pivot_wider(names_from = cohort, values_from = n, values_fill = 0)

totals <- long_counts %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    total_microbiome = sum(.data$autoimmune_microbiome, na.rm = TRUE),
    total_diet       = sum(.data$autoimmune_diet,       na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- long_counts %>%
  dplyr::left_join(totals, by = "variable") %>%
  tidyr::pivot_longer(
    cols = c(autoimmune_microbiome, autoimmune_diet),
    names_to = "cohort", values_to = "count"
  ) %>%
  dplyr::mutate(
    total_cohort = dplyr::case_when(
      cohort == "autoimmune_microbiome" ~ total_microbiome,
      cohort == "autoimmune_diet"       ~ total_diet,
      TRUE ~ NA_real_
    ),
    prop = dplyr::if_else(total_cohort > 0, count / total_cohort, 0)
  )

# ---------- 5) plot (namespaced ggplot2 calls; no library(ggplot2) needed) ----------
p <- ggplot2::ggplot(plot_df_top, ggplot2::aes(x = value_plot, y = prop, fill = cohort)) +
  ggplot2::geom_col(width = 0.9, color = "gray30", linewidth = 0.15) +
  ggplot2::facet_wrap(~ variable, scales = "free_x", ncol = 2) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                              expand = ggplot2::expansion(mult = c(0, 0.02))) +
  ggplot2::scale_fill_brewer(palette = "Set2", name = "Cohort") +
  ggplot2::labs(
    title = "Key SRA Classifications for Biomarker Robustness",
    subtitle = "Proportion within each variable; top categories (others lumped)",
    x = "Category",
    y = "Proportion within variable"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
    axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text         = ggplot2::element_text(face = "bold")
  )

print(p)
ggplot2::ggsave("sra_biomarker_key_classifications_stacked_100pct.png", p, width = 12, height = 9, dpi = 300)
ggplot2::ggsave("sra_biomarker_key_classifications_stacked_100pct.pdf", p, width = 12, height = 9)
``
