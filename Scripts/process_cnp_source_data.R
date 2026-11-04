# =============================================================================
# CuratedNutritionProfiler (CNP) — process downloaded GEO data into compact,
# analysis-ready tables.
#
# Input : data-raw/CNP_Source_Data/<group>/<accession>/...   (from
#         Scripts/download_cnp_source_data.R)
# Output: data/processed/<group>/<dataset>_expr.csv.gz
#             gene-symbol x sample matrix on a log2 scale
#             (log2 CPM for counts; log2(x + 1) for FPKM/TPM; log2 for
#             arrays reported on a linear scale; arrays already in log2 kept)
#         data/processed/<group>/<dataset>_pheno.csv
#             one row per sample: sample, condition (original study arm),
#             plus study covariates (genotype, batch, time, ...)
#         data/processed/dataset_index.csv
#             one row per processed dataset (species, platform, units, n, ...)
#
# Case/control labels are NOT assigned here. Which conditions count as
# "low nutrient" vs "replete" is a modelling decision, defined in one place in
# Reports/CNP_Signature_Validation.Rmd (object `panels`).
#
# Gene identifiers are converted to the dataset's native gene symbols
# (org.Hs.eg.db / org.Mm.eg.db, or the GEO platform annotation for arrays).
# Several probes/IDs per symbol are collapsed by keeping the one with the
# highest mean expression.
#
# Not processed (no usable nutrient contrast or not expression data):
#   GSE266352 (LAMN vs CRC metastases), GSE124395 (liver cell atlas),
#   GSE246328 (NASH diet, no iron arm), GSE296260 (ccRCC vs kidney),
#   GSE143307 (DNA methylation), GSE271214/GSE271249 (2 RA-treated samples
#   per arm, confounded with genotype/day; RAW tar not downloaded),
#   GDS4860 (curated subset of GSE46268, which is processed instead).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Matrix)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
})

PROJECT_DIR <- if (file.exists("curatedNutritionDataset.Rproj")) getwd() else
  path.expand("~/RStudio Projects/curatedNutritionDataset")
RAW  <- file.path(PROJECT_DIR, "data-raw", "CNP_Source_Data")
OUT  <- file.path(PROJECT_DIR, "data", "processed")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----------------------------------------------------------------
raw_path <- function(...) file.path(RAW, ...)

read_series_matrix <- function(file) {
  lines <- readLines(gzfile(file))
  meta  <- grep("^!Sample_", lines, value = TRUE)
  split_line <- function(l) gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]])
  m <- lapply(meta, split_line)
  keys <- vapply(m, `[`, "", 1)
  vals <- lapply(m, `[`, -1)
  ph <- data.frame(gsm = vals[[which(keys == "!Sample_geo_accession")[1]]],
                   title = vals[[which(keys == "!Sample_title")[1]]],
                   source = vals[[which(keys == "!Sample_source_name_ch1")[1]]],
                   stringsAsFactors = FALSE)
  for (j in which(keys == "!Sample_characteristics_ch1")) {
    v <- vals[[j]]; k <- sub(":.*", "", v[1])
    ph[[make.names(k)]] <- trimws(sub("^[^:]*:", "", v))
  }
  b <- grep("^!series_matrix_table_begin", lines); e <- grep("^!series_matrix_table_end", lines)
  expr <- if (length(b) && e - b > 2) {
    x <- fread(text = lines[(b + 1):(e - 1)], sep = "\t", quote = "\"", data.table = FALSE)
    rn <- as.character(x[[1]]); x <- as.matrix(x[, -1, drop = FALSE]); rownames(x) <- rn
    storage.mode(x) <- "double"; x
  } else NULL
  list(pheno = ph, expr = expr)
}

read_soft_platform <- function(soft_file) {
  con <- gzfile(soft_file); on.exit(close(con))
  lines <- readLines(con)
  b <- grep("^!platform_table_begin", lines)[1]; e <- grep("^!platform_table_end", lines)[1]
  fread(text = lines[(b + 1):(e - 1)], sep = "\t", quote = "", data.table = FALSE, fill = TRUE)
}

symbol_col <- function(pl) {
  for (k in c("Gene Symbol", "GENE_SYMBOL", "Symbol", "ORF")) if (k %in% names(pl)) return(k)
  stop("no symbol column in platform table")
}

# Keep, for each symbol, the row with the highest mean expression
collapse_to_symbol <- function(mat, symbols) {
  symbols <- trimws(sub(" ///.*", "", as.character(symbols)))
  ok <- !is.na(symbols) & symbols != "" & symbols != "---"
  mat <- mat[ok, , drop = FALSE]; symbols <- symbols[ok]
  o <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[o, , drop = FALSE]; symbols <- symbols[o]
  keep <- !duplicated(symbols)
  out <- mat[keep, , drop = FALSE]; rownames(out) <- symbols[keep]
  out[order(rownames(out)), , drop = FALSE]
}

map_ids <- function(ids, species, keytype) {
  db <- if (species == "Homo sapiens") org.Hs.eg.db else org.Mm.eg.db
  suppressMessages(mapIds(db, keys = as.character(ids), column = "SYMBOL",
                          keytype = keytype, multiVals = "first"))
}

# Library-size QC for count data: drop samples whose total counts are < 10% of the
# median library (failed libraries); returns log2 CPM of the remaining samples
qc_dropped <- character(0)
log_cpm <- function(counts) {
  lib <- colSums(counts, na.rm = TRUE)
  bad <- lib < 0.1 * median(lib)
  if (any(bad)) {
    message(sprintf("  QC: dropping %s (library %s reads vs median %s)", paste(colnames(counts)[bad], collapse = ", "),
                    paste(format(lib[bad], big.mark = ","), collapse = ", "), format(median(lib), big.mark = ",")))
    qc_dropped <<- c(qc_dropped, colnames(counts)[bad])
    counts <- counts[, !bad, drop = FALSE]
  }
  edgeR::cpm(counts, log = TRUE, prior.count = 1)
}
log_plus <- function(x) log2(x + 1)
log_array <- function(x) if (max(x, na.rm = TRUE) > 100) log2(pmax(x, 1)) else x

# Run a dataset block only when its raw download is present. Raw files are deleted after
# processing, so a re-run keeps the existing processed output (and index row) instead.
skipped_acc <- character(0)
with_raw <- function(group, acc, expr) {
  if (!dir.exists(raw_path(group, acc))) {
    message(sprintf("  %-20s raw files not present - keeping existing processed output", acc))
    skipped_acc <<- c(skipped_acc, acc); return(invisible())
  }
  eval(substitute(expr), envir = new.env(parent = globalenv()))
}

# Gene symbol from Affymetrix "gene_assignment"-style strings: "NM_x // SYMBOL // name // ..."
parse_assignment <- function(x) {
  s <- trimws(sub(" *//.*$", "", sub("^[^/]*// *", "", as.character(x))))
  s[is.na(x) | s %in% c("", "---")] <- NA; s
}
# Clontech Atlas arrays (e.g. GPL155) give only names such as "cyclin G-associated kinase (GAK)":
# take the last parenthesised symbol and write it in rodent case (Gak)
parse_paren_symbol <- function(x) {
  m <- regmatches(x, gregexpr("\\(([A-Za-z0-9-]+)\\)", x))
  s <- vapply(m, function(v) if (length(v)) gsub("[()]", "", tail(v, 1)) else NA_character_, "")
  ifelse(is.na(s), NA, paste0(substr(s, 1, 1), tolower(substring(s, 2))))
}

index <- list()
save_dataset <- function(id, accession, group, species, platform, units, expr, pheno, note = "") {
  # samples removed by QC are dropped from the phenotype table and recorded in the note
  dropped <- setdiff(pheno$sample, colnames(expr))
  if (length(dropped)) note <- paste(c(note[nzchar(note)], paste("QC-excluded (failed library):",
                                       paste(dropped, collapse = ", "))), collapse = "; ")
  pheno <- pheno[pheno$sample %in% colnames(expr), , drop = FALSE]
  expr  <- expr[, pheno$sample, drop = FALSE]
  stopifnot(identical(colnames(expr), pheno$sample))
  expr <- expr[rowSums(is.finite(expr)) == ncol(expr), , drop = FALSE]
  d <- file.path(OUT, group); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ef <- file.path(d, paste0(id, "_expr.csv.gz")); pf <- file.path(d, paste0(id, "_pheno.csv"))
  dt <- data.table(gene = rownames(expr), round(expr, 4))
  setnames(dt, c("gene", colnames(expr)))
  fwrite(dt, ef, compress = "gzip")
  fwrite(pheno, pf)
  index[[id]] <<- data.frame(dataset = id, accession = accession, group = group, species = species,
                             platform = platform, units = units, n_samples = ncol(expr),
                             n_genes = nrow(expr),
                             conditions = paste(sprintf("%s (%d)", names(table(pheno$condition)),
                                                        as.vector(table(pheno$condition))), collapse = "; "),
                             expr_file = sub(paste0(PROJECT_DIR, "/"), "", ef, fixed = TRUE),
                             pheno_file = sub(paste0(PROJECT_DIR, "/"), "", pf, fixed = TRUE),
                             note = note, stringsAsFactors = FALSE)
  message(sprintf("  %-20s %5d genes x %3d samples  -> %s", id, nrow(expr), ncol(expr),
                  format(structure(file.size(ef), class = "object_size"), units = "auto")))
}

array_dataset <- function(group, acc, symbol_from = c("platform", "entrez"), sym_fun = NULL,
                          force_log2 = FALSE) {
  symbol_from <- match.arg(symbol_from)
  sm <- read_series_matrix(list.files(raw_path(group, acc, "matrix"), full.names = TRUE)[1])
  x  <- if (force_log2) log2(pmax(sm$expr, 1e-3)) else log_array(sm$expr)
  sym <- if (symbol_from == "platform") {
    pl <- read_soft_platform(raw_path(group, acc, "soft", paste0(acc, "_family.soft.gz")))
    sy <- if (is.null(sym_fun)) pl[[symbol_col(pl)]] else sym_fun(pl)
    sy[match(rownames(x), as.character(pl$ID))]
  } else {
    map_ids(sub("_at$", "", rownames(x)), "Homo sapiens", "ENTREZID")
  }
  list(expr = collapse_to_symbol(x, sym), pheno = sm$pheno)
}

ncbi_counts_dataset <- function(group, acc) {   # human series only (GRCh38)
  f <- raw_path(group, acc, "ncbi_counts", paste0(acc, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
  x <- fread(f, data.table = FALSE); ids <- x$GeneID
  x <- as.matrix(x[, -1]); rownames(x) <- ids
  lc <- log_cpm(x)
  list(expr = collapse_to_symbol(lc, map_ids(ids, "Homo sapiens", "ENTREZID")),
       pheno = read_series_matrix(list.files(raw_path(group, acc, "matrix"), full.names = TRUE)[1])$pheno)
}

message("==== Processing CNP source data -> ", OUT)

# ---- Vitamin A --------------------------------------------------------------
# GSE134933: mouse kidney (WT / TRACK x chow / VAD) and liver (WT chow / VAD), FPKM
with_raw("Vitamin_A", "GSE134933", {
  tf <- tempfile(fileext = ".xls")
  R.utils::gunzip(raw_path("Vitamin_A", "GSE134933", "suppl",
                           "GSE134933_2018-2-27_expression.genes.protein_coding.xls.gz"),
                  tf, remove = FALSE)
  x <- as.data.frame(read_excel(tf, sheet = "gene"))
  smp <- setdiff(names(x), c("Tracking_id", "Gene_id", "Gene_name", "Locus", "biotype", "ensembl_id"))
  m <- as.matrix(x[, smp]); storage.mode(m) <- "double"
  m <- collapse_to_symbol(log_plus(m), x$Gene_name)
  kid <- grep("^liver", smp, invert = TRUE, value = TRUE)
  liv <- grep("^liver", smp, value = TRUE)
  ph_k <- data.frame(sample = kid,
                     condition = ifelse(grepl("VAD", kid), "VAD", "Chow"),
                     genotype  = ifelse(grepl("TRACK", kid), "TRACK", "WT"))
  save_dataset("GSE134933_kidney", "GSE134933", "Vitamin_A", "Mus musculus",
               "RNA-seq (GPL13112)", "log2(FPKM+1)", m[, kid], ph_k)
  ph_l <- data.frame(sample = liv, condition = ifelse(grepl("vad", liv), "VAD", "Chow"))
  save_dataset("GSE134933_liver", "GSE134933", "Vitamin_A", "Mus musculus",
               "RNA-seq (GPL13112)", "log2(FPKM+1)", m[, liv], ph_l,
               "Liver samples deposited in the same series as the kidney data")
})

# GSE299530: mouse lung VAD vs VAS, counts keyed by Entrez ID
with_raw("Vitamin_A", "GSE299530", {
  x <- fread(raw_path("Vitamin_A", "GSE299530", "suppl", "GSE299530_features_combined_2021.txt.gz"),
             data.table = FALSE)
  ids <- x[[1]]; x <- as.matrix(x[, -1]); rownames(x) <- ids
  # Columns from 2019 are the BALF samples (batch 1, mostly microbial reads) - lung only here
  x <- x[, !grepl("_2019$", colnames(x))]
  m <- collapse_to_symbol(log_cpm(x), map_ids(ids, "Mus musculus", "ENTREZID"))
  s <- colnames(m)
  ph <- data.frame(sample = s, condition = sub("\\d+_.*", "", s),
                   timepoint = sub("^[^_]+_([^_]+)_.*", "\\1", s),
                   batch = sub(".*_", "", s))
  save_dataset("GSE299530", "GSE299530", "Vitamin_A", "Mus musculus",
               "RNA-seq (GPL35938/GPL35939)", "log2 CPM", m, ph,
               "BALF samples (2019 batch) excluded; whole-lung only")
})

# GSE37313: rat liver VAD/VAS x SHAM/BDL x vehicle/retinyl palmitate, Illumina array
with_raw("Vitamin_A", "GSE37313", {
  a <- array_dataset("Vitamin_A", "GSE37313")
  p <- a$pheno
  ph <- data.frame(sample = p$gsm,
                   condition = paste(ifelse(grepl("^VAD|deficient", p$treatment.1), "VAD", "VAS"),
                                     ifelse(grepl("VITA", p$treatment.3), "RetPalm", "Vehicle"), sep = "_"),
                   surgery = ifelse(grepl("BDL", p$treatment.2), "BDL", "SHAM"))
  save_dataset("GSE37313", "GSE37313", "Vitamin_A", "Rattus norvegicus",
               "Illumina RatRef-12 (GPL6101)", "log2 intensity", a$expr[, ph$sample], ph)
})

# GSE46268: human monocytes, control vs ATRA vs 1,25D3 (4 donors), Affymetrix GPL570
with_raw("Vitamin_A", "GSE46268", {
  a <- array_dataset("Vitamin_A", "GSE46268")
  p <- a$pheno
  cond <- c("control" = "Control", "all-trans retinoic acid" = "ATRA",
            "1,25a-dihydroxyvitamin D3" = "Calcitriol")[p$treatment]
  ph <- data.frame(sample = p$gsm, condition = unname(cond), donor = p$individual)
  save_dataset("GSE46268", "GSE46268", "Vitamin_A", "Homo sapiens",
               "Affymetrix HG-U133 Plus 2 (GPL570)", "log2 MAS5", a$expr[, ph$sample], ph,
               "Used for both the vitamin A (ATRA) and vitamin D (calcitriol) panels")
})

# ---- Vitamin D --------------------------------------------------------------
with_raw("Vitamin_D", "GSE63086", {   # GSE63086: MG63 cells, vehicle vs 1,25D x control vs HNRNPC knock-down
  d <- ncbi_counts_dataset("Vitamin_D", "GSE63086"); p <- d$pheno
  ph <- data.frame(sample = p$gsm,
                   condition = ifelse(grepl("^Vitamin D", p$title), "Calcitriol", "Vehicle"),
                   knockdown = ifelse(grepl("KD", p$title), "HNRNPC_KD", "Control"))
  save_dataset("GSE63086", "GSE63086", "Vitamin_D", "Homo sapiens", "RNA-seq (GPL16791)",
               "log2 CPM (NCBI counts)", d$expr[, intersect(ph$sample, colnames(d$expr))], ph)
})

with_raw("Vitamin_D", "GSE221690", {   # GSE221690: 4 prostate cell lines, vehicle vs 1,25D3
  d <- ncbi_counts_dataset("Vitamin_D", "GSE221690"); p <- d$pheno
  ph <- data.frame(sample = p$gsm,
                   condition = ifelse(grepl("treat", p$title), "Calcitriol", "Vehicle"),
                   cell_line = p$cell.line)
  save_dataset("GSE221690", "GSE221690", "Vitamin_D", "Homo sapiens", "RNA-seq (GPL16791)",
               "log2 CPM (NCBI counts)", d$expr[, intersect(ph$sample, colnames(d$expr))], ph)
})

with_raw("Vitamin_D", "GSE15947", {   # GSE15947: RWPE1 cells, vehicle vs 1,25D at 6/24/48 h, GPL570 (already log2)
  a <- array_dataset("Vitamin_D", "GSE15947"); p <- a$pheno
  ph <- data.frame(sample = p$gsm,
                   condition = ifelse(grepl("Vehicle", p$treatment), "Vehicle", "Calcitriol"),
                   time = p$time.point)
  save_dataset("GSE15947", "GSE15947", "Vitamin_D", "Homo sapiens",
               "Affymetrix HG-U133 Plus 2 (GPL570)", "log2 intensity (RMA)", a$expr[, ph$sample], ph)
})

with_raw("Vitamin_D", "GSE283020", {   # GSE283020: mouse TH2 cells, regular vs vitamin D-deficient diet, counts keyed by Ensembl
  x <- fread(raw_path("Vitamin_D", "GSE283020", "suppl", "GSE283020_mouse_vitd_def_th2_matrix.csv.gz"),
             data.table = FALSE)
  ids <- x[[1]]; x <- as.matrix(x[, -1]); rownames(x) <- ids
  m <- collapse_to_symbol(log_cpm(x), map_ids(ids, "Mus musculus", "ENSEMBL"))
  p <- read_series_matrix(list.files(raw_path("Vitamin_D", "GSE283020", "matrix"), full.names = TRUE)[1])$pheno
  p <- p[match(colnames(m), p$title), ]
  ph <- data.frame(sample = colnames(m),
                   condition = ifelse(grepl("deficient", p$treatment), "VitD_deficient", "Regular"),
                   batch = p$batch)
  save_dataset("GSE283020", "GSE283020", "Vitamin_D", "Mus musculus", "RNA-seq (GPL17021)",
               "log2 CPM", m, ph)
})

# ---- Vitamin C --------------------------------------------------------------
with_raw("Vitamin_C", "GSE233598", {   # GSE233598: mouse liver, Gulo-/- vitamin C dose series + WT, TPM
  x <- as.data.frame(read_excel(raw_path("Vitamin_C", "GSE233598", "suppl", "GSE233598_tpm.xlsx"), sheet = "tpm"))
  smp <- setdiff(names(x), c("ENSEMBL", "SYMBOL"))
  m <- as.matrix(x[, smp]); storage.mode(m) <- "double"
  m <- collapse_to_symbol(log_plus(m), x$SYMBOL)
  code <- sub("^[FM]", "", sub("-.*", "", smp))   # GL00, GL01, GL05, GL40, GLR40, WT00
  lab <- c(GL00 = "Gulo_0pct", GL01 = "Gulo_0.01pct", GL05 = "Gulo_0.05pct",
           GL40 = "Gulo_0.40pct", GLR40 = "Gulo_rescue_0.40pct", WT00 = "WT")
  ph <- data.frame(sample = smp, condition = unname(lab[code]),
                   sex = ifelse(startsWith(smp, "F"), "Female", "Male"))
  save_dataset("GSE233598", "GSE233598", "Vitamin_C", "Mus musculus", "RNA-seq (GPL24247)",
               "log2(TPM+1)", m, ph)
})

with_raw("Vitamin_C", "GSE190371", {   # GSE190371: mouse kidney scRNA-seq, Gulo-/- control vs vitamin C-deficient -> pseudobulk
  tarf <- raw_path("Vitamin_C", "GSE190371", "suppl", "GSE190371_RAW.tar")
  td <- tempfile(); dir.create(td)
  members <- grep("(barcodes|genes|matrix)\\.(tsv|mtx)\\.gz$", untar(tarf, list = TRUE), value = TRUE)
  untar(tarf, files = members, exdir = td)
  stems <- unique(sub("_(barcodes|genes|matrix)\\..*$", "", basename(members)))
  pb <- NULL; gene_ids <- NULL
  for (s in stems) {
    mm <- readMM(file.path(td, paste0(s, "_matrix.mtx.gz")))
    g  <- fread(file.path(td, paste0(s, "_genes.tsv.gz")), header = FALSE, data.table = FALSE)
    v  <- Matrix::rowSums(mm); names(v) <- g[[1]]
    if (is.null(gene_ids)) gene_ids <- g[, 1:2]
    pb <- cbind(pb, v[gene_ids[[1]]])
  }
  colnames(pb) <- sub("^GSM\\d+_", "", stems); rownames(pb) <- gene_ids[[1]]
  unlink(td, recursive = TRUE)
  m <- collapse_to_symbol(log_cpm(pb), gene_ids[[2]])
  ph <- data.frame(sample = colnames(m),
                   condition = ifelse(grepl("^vcd", colnames(m), ignore.case = TRUE), "VitC_deficient", "Control"))
  save_dataset("GSE190371", "GSE190371", "Vitamin_C", "Mus musculus",
               "scRNA-seq pseudobulk (GPL21273)", "log2 CPM (summed UMIs)", m, ph,
               "Pseudobulk of 2 control vs 2 deficient kidneys - descriptive only")
})

# ---- Vitamin B12 ------------------------------------------------------------
with_raw("Vitamin_B12", "GSE161763", {   # GSE161763: mouse embryonic brain, Ronin F80L (cblX) vs control, counts
  x <- fread(raw_path("Vitamin_B12", "GSE161763", "suppl", "GSE161763_RONIN_F80L_RNAseq_raw_counts.txt.gz"),
             data.table = FALSE)
  g <- x$Gene; x <- as.matrix(x[, -1]); rownames(x) <- g
  m <- collapse_to_symbol(log_cpm(x), g)
  ph <- data.frame(sample = colnames(m), condition = ifelse(grepl("mutant", colnames(m)), "cblX_F80L", "Control"))
  save_dataset("GSE161763", "GSE161763", "Vitamin_B12", "Mus musculus", "RNA-seq (GPL17021)", "log2 CPM", m, ph)
})

with_raw("Vitamin_B12", "GSE103417", {   # GSE103417: N1E-115 cells, WT / OT (no Cbl binding) / TO (sequesters Cbl), Agilent
  a <- array_dataset("Vitamin_B12", "GSE103417"); p <- a$pheno
  arm <- sub("^(proliferated|differentiated)_neuroblastoma_", "", p$source)
  ph <- data.frame(sample = p$gsm,
                   condition = c(WT = "WT", OT = "OT_control_fusion", TO = "TO_B12_sequestered")[arm],
                   medium = p$medium)
  save_dataset("GSE103417", "GSE103417", "Vitamin_B12", "Mus musculus",
               "Agilent Whole Mouse Genome (GPL4134)", "log2 intensity", a$expr[, ph$sample], ph)
})

with_raw("Vitamin_B12", "GSE104164", {   # GSE104164: rat cerebellum, control vs methyl-donor (B12/folate) deficiency, Agilent
  a <- array_dataset("Vitamin_B12", "GSE104164"); p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("^MDD", p$source), "MDD", "Control"))
  save_dataset("GSE104164", "GSE104164", "Vitamin_B12", "Rattus norvegicus",
               "Agilent Rat Gene Expression (GPL14797)", "log2 intensity", a$expr[, ph$sample], ph)
})

with_raw("Vitamin_B12", "GSE99113", {   # GSE99113: mouse astrocytes, B12-sufficient vs B12-free medium (2 vs 2)
  x <- fread(raw_path("Vitamin_B12", "GSE99113", "suppl", "GSE99113_JC_NS190_AstrocyteB12.txt.gz"),
             data.table = FALSE)
  g <- x$name; x <- as.matrix(x[, -1]); rownames(x) <- g
  m <- collapse_to_symbol(log_plus(x), g)
  ph <- data.frame(sample = colnames(m), condition = ifelse(grepl("def", colnames(m)), "B12_deficient", "B12_sufficient"))
  save_dataset("GSE99113", "GSE99113", "Vitamin_B12", "Mus musculus", "RNA-seq (GPL19057)",
               "log2(normalised expression+1)", m, ph, "2 vs 2 - descriptive only")
})

# ---- Vitamin K --------------------------------------------------------------
with_raw("Vitamin_K", "GSE148011", {   # GSE148011: p53/PTEN-null MEFs, vitamin K precursor (MSB) 0/10/20 uM (+ erastin arms)
  tarf <- raw_path("Vitamin_K", "GSE148011", "suppl", "GSE148011_RAW.tar")
  td <- tempfile(); dir.create(td); untar(tarf, exdir = td)
  p <- read_series_matrix(list.files(raw_path("Vitamin_K", "GSE148011", "matrix"), full.names = TRUE)[1])$pheno
  per_sample <- lapply(p$gsm, function(gsm) {
    f <- list.files(td, pattern = paste0("^", gsm, "_"), full.names = TRUE)
    y <- fread(f, header = FALSE, data.table = FALSE)
    y <- y[!grepl("^__", y[[1]]), ]
    tapply(y[[2]], sub("\\.\\d+$", "", y[[1]]), sum)      # strip Ensembl version
  })
  # The vitamin K and erastin arms were quantified against different gene annotations
  # (45,323 vs 53,384 IDs): keep only genes present in every file rather than zero-filling
  genes <- sort(Reduce(intersect, lapply(per_sample, names)))
  cnt <- vapply(per_sample, function(v) v[genes], numeric(length(genes)))
  rownames(cnt) <- genes
  unlink(td, recursive = TRUE); colnames(cnt) <- p$gsm
  m <- collapse_to_symbol(log_cpm(cnt), map_ids(rownames(cnt), "Mus musculus", "ENSEMBL"))
  ph <- data.frame(sample = p$gsm, condition = sub("_n\\d+$", "", p$title))   # K3_0, K3_10, K3_20, Erastin_0, ...
  save_dataset("GSE148011", "GSE148011", "Vitamin_K", "Mus musculus", "RNA-seq (GPL19057)", "log2 CPM", m, ph)
})

# ---- Omega-3 ----------------------------------------------------------------
with_raw("Omega_3", "GSE12375", {   # GSE12375: human PBMC, EPA/DHA vs MUFA (oleic acid) at baseline and 6 months, NuGO array
  a <- array_dataset("Omega_3", "GSE12375", symbol_from = "entrez"); p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = sub(".*_(EPADHA|MUFA)_(BASE|END)$", "\\1_\\2", p$title),
                   subject = sub("^A101_([^_]+)_.*", "\\1", p$title))
  save_dataset("GSE12375", "GSE12375", "Omega_3", "Homo sapiens",
               "NuGO Hs1a520180 custom Affymetrix (GPL7144)", "log2 intensity", a$expr[, ph$sample], ph)
})

# ---- Vitamin B6 (added 2026-11-04) -------------------------------------------
with_raw("Vitamin_B6", "GSE345689", {   # mouse kidney, low vs high B6 diet x WT / Sgpl1 R222Q; FPKM in DE workbooks
  d <- raw_path("Vitamin_B6", "GSE345689", "suppl")
  rd <- function(f) { x <- as.data.frame(read_excel(file.path(d, f)));
    x[, c("gene_id", "gene_name", grep("^FPKM\\.", names(x), value = TRUE))] }
  a <- rd("GSE345689_WT_HVSWT_L_Gene_differential_expression.xlsx")       # WT high + WT low
  b <- rd("GSE345689_KI_HVSKI_L_Gene_differential_expression.xlsx")       # KI high + KI low
  x <- merge(a, b[, c("gene_id", grep("^FPKM\\.", names(b), value = TRUE))], by = "gene_id")
  m <- as.matrix(x[, grep("^FPKM\\.", names(x))]); storage.mode(m) <- "double"
  colnames(m) <- sub("^FPKM\\.", "", colnames(m))
  m <- collapse_to_symbol(log_plus(m), x$gene_name)
  ph <- data.frame(sample = colnames(m),
                   condition = ifelse(grepl("low", colnames(m), ignore.case = TRUE), "Low_B6", "High_B6"),
                   genotype = ifelse(startsWith(colnames(m), "KI"), "Sgpl1_R222Q", "WT"))
  save_dataset("GSE345689", "GSE345689", "Vitamin_B6", "Mus musculus", "RNA-seq (GPL24247)",
               "log2(FPKM+1)", m, ph, "Expression taken from the deposited DE workbooks (the 'raw_counts' file holds read QC only)")
})

with_raw("Vitamin_B6", "GSE139892", {   # mouse leukemic cells, PDXK inhibitor 4-O-methylpyridoxine vs vehicle
  x <- fread(raw_path("Vitamin_B6", "GSE139892", "suppl", "GSE139892_Chen_RNASeq_OMP_RPKM.csv.gz"), data.table = FALSE)
  m <- as.matrix(x[, grep("Day6", names(x))]); m <- collapse_to_symbol(log_plus(m), x$Gene_name)
  ph <- data.frame(sample = colnames(m), condition = ifelse(grepl("^OMP", colnames(m)), "OMP_B6_antagonist", "Vehicle"))
  save_dataset("GSE139892", "GSE139892", "Vitamin_B6", "Mus musculus", "RNA-seq", "log2(RPKM+1)", m, ph)
})

# ---- Folate (added 2026-11-04) ----------------------------------------------
with_raw("Folate", "GSE275393", {   # E11.5 brain, maternal folic acid 0.1 / 3 / 30 ppm x Tet1 genotype, counts
  x <- fread(raw_path("Folate", "GSE275393", "suppl", "GSE275393_E11-5_brains_FA_counts.txt.gz"), data.table = FALSE)
  smp <- names(x)[-(1:4)]; m <- as.matrix(x[, smp]); rownames(m) <- x$ensembl_gene_id
  m <- collapse_to_symbol(log_cpm(m), x$mgi_symbol)
  ppm <- sub(".*_(01|3|30)ppm_.*", "\\1", colnames(m))
  ph <- data.frame(sample = colnames(m),
                   condition = c("01" = "FA_0.1ppm", "3" = "FA_3ppm", "30" = "FA_30ppm")[ppm],
                   genotype = sub(".*ppm_(WT|HET|KO)_.*", "Tet1_\\1", colnames(m)),
                   sex = sub(".*ppm_(WT|HET|KO)_([FM])_.*", "\\2", colnames(m)))
  save_dataset("GSE275393", "GSE275393", "Folate", "Mus musculus", "RNA-seq", "log2 CPM", m, ph)
})

with_raw("Folate", "GSE324232", {   # pachytene spermatocytes + round spermatids, folate-deficient vs sufficient, counts
  rd <- function(f) { x <- fread(raw_path("Folate", "GSE324232", "suppl", f), data.table = FALSE)
    m <- as.matrix(x[, -1]); rownames(m) <- sub("\\.\\d+$", "", x[[1]]); m }
  a <- rd("GSE324232_WT_PS_RNA_counts_matrix.csv.gz"); b <- rd("GSE324232_WT_RS_RNA_counts_matrix.csv.gz")
  g <- intersect(rownames(a), rownames(b)); cnt <- cbind(a[g, ], b[g, ])
  m <- collapse_to_symbol(log_cpm(cnt), map_ids(rownames(cnt), "Mus musculus", "ENSEMBL"))
  ph <- data.frame(sample = colnames(m),
                   condition = ifelse(grepl("FD", colnames(m)), "Folate_deficient", "Folate_sufficient"),
                   cell_type = ifelse(startsWith(colnames(m), "PS"), "Pachytene_spermatocyte", "Round_spermatid"))
  save_dataset("GSE324232", "GSE324232", "Folate", "Mus musculus", "RNA-seq", "log2 CPM", m, ph,
               "Labels from count-matrix column names (3 FD + 3 FS per cell type); GEO characteristics list 5 FD / 7 FS")
})

with_raw("Folate", "GSE271941", {   # MEFs, moderate vs fortified folic acid x Rsg1 genotype; columns follow GEO sample order
  x <- fread(raw_path("Folate", "GSE271941", "suppl", "GSE271941_salmon.merged.gene_counts.tsv.gz"), data.table = FALSE)
  m <- as.matrix(x[, -(1:2)]); rownames(m) <- x$gene_id
  m <- collapse_to_symbol(log_cpm(round(m)), x$gene_name)
  # GEO: "Processed matrix column 1-6" = moderate FA (BHOM, BWT), 7-12 = fortified FA (PHOM, PWT)
  ph <- data.frame(sample = colnames(m),
                   condition = ifelse(startsWith(colnames(m), "B"), "Moderate_FA", "Fortified_FA"),
                   genotype = ifelse(grepl("HOM", colnames(m)), "Rsg1_L3P", "WT"))
  save_dataset("GSE271941", "GSE271941", "Folate", "Mus musculus", "RNA-seq", "log2 CPM", m, ph,
               "Folic acid dose contrast (moderate vs fortified), not deficiency")
})

# ---- Carnitine (added 2026-11-04) -------------------------------------------
with_raw("Carnitine", "GSE343742", {   # BMM +/- L-carnitine, +/- LPS priming, TPM
  x <- fread(raw_path("Carnitine", "GSE343742", "suppl", "GSE343742_20250115_Elliott_RNAseq_tpm.csv.gz"), data.table = FALSE)
  smp <- grep("^Neg", names(x), value = TRUE)
  m <- collapse_to_symbol(log_plus(as.matrix(x[, smp])), x$symbol)
  arm <- sub("_\\d+$", "", colnames(m))
  ph <- data.frame(sample = colnames(m),
                   condition = c(NegCtrl = "Untreated", NegLC6 = "Carnitine", NegLPS = "LPS", NegLPSLC6 = "LPS_Carnitine")[arm],
                   lps = ifelse(grepl("LPS", arm), "LPS", "None"))
  save_dataset("GSE343742", "GSE343742", "Carnitine", "Mus musculus", "RNA-seq", "log2(TPM+1)", m, ph)
})

with_raw("Carnitine", "GSE344237", {   # liver, control vs high-fat diet x resveratrol + L-carnitine (none/low/high), Agilent
  a <- array_dataset("Carnitine", "GSE344237"); p <- a$pheno
  arm <- trimws(sub(".* - ", "", p$title))
  ph <- data.frame(sample = p$gsm, condition = arm,
                   diet = ifelse(startsWith(arm, "HFCD"), "HFCD", "CD"),
                   supplement = ifelse(grepl("RCh", arm), "High", ifelse(grepl("RCl", arm), "Low", "None")))
  save_dataset("GSE344237", "GSE344237", "Carnitine", "Mus musculus", "Agilent SurePrint G3 Mouse (GPL21163)",
               "log2 intensity", a$expr[, intersect(ph$sample, colnames(a$expr))], ph,
               "L-carnitine always co-administered with resveratrol")
})

# ---- Vitamin K external sets (added 2026-11-04) ------------------------------
with_raw("Vitamin_K", "GSE267133", {   # ARPE-19, erastin vs erastin + vitamin K1, per-sample FPKM files
  tarf <- raw_path("Vitamin_K", "GSE267133", "suppl", "GSE267133_RAW.tar")
  td <- tempfile(); dir.create(td); untar(tarf, exdir = td)
  p <- read_series_matrix(list.files(raw_path("Vitamin_K", "GSE267133", "matrix"), full.names = TRUE)[1])$pheno
  vals <- lapply(p$gsm, function(gsm) {
    y <- fread(list.files(td, pattern = paste0("^", gsm, "_"), full.names = TRUE), data.table = FALSE)
    setNames(y[[grep("^FPKM", names(y))]], y$gene_name) })
  unlink(td, recursive = TRUE)
  g <- Reduce(intersect, lapply(vals, names))
  m <- vapply(vals, function(v) tapply(v, names(v), max)[g], numeric(length(g))); colnames(m) <- p$gsm
  m <- collapse_to_symbol(log_plus(m), g)
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("vk", p$title), "Erastin_VitK1", "Erastin"))
  save_dataset("GSE267133", "GSE267133", "Vitamin_K", "Homo sapiens", "RNA-seq", "log2(FPKM+1)", m, ph)
})

with_raw("Vitamin_K", "GSE237219", {   # Jurkat, vitamin K1 vs control (NCBI counts)
  d <- ncbi_counts_dataset("Vitamin_K", "GSE237219"); p <- d$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("^VK1", p$title), "VitK1", "Control"))
  save_dataset("GSE237219", "GSE237219", "Vitamin_K", "Homo sapiens", "RNA-seq",
               "log2 CPM (NCBI counts)", d$expr[, intersect(ph$sample, colnames(d$expr))], ph)
})

with_raw("Vitamin_K", "GSE272401", {   # mouse BMDM: PBS / RANKL / RANKL + vitamins D and K, FPKM
  x <- fread(raw_path("Vitamin_K", "GSE272401", "suppl", "GSE272401_genes_fpkm_expression_updated_240717.txt.gz"),
             data.table = FALSE, quote = "")
  smp <- grep("^FPKM\\.", names(x), value = TRUE)
  m <- as.matrix(x[, smp]); colnames(m) <- sub("^FPKM\\.", "", smp)
  m <- collapse_to_symbol(log_plus(m), x$gene_name)
  arm <- sub("_\\d+$", "", colnames(m))
  ph <- data.frame(sample = colnames(m),
                   condition = c(Rankl = "RANKL", Negative = "PBS", Rankl_Vit_D_K = "RANKL_VitD_VitK")[arm])
  save_dataset("GSE272401", "GSE272401", "Vitamin_K", "Mus musculus", "RNA-seq", "log2(FPKM+1)", m, ph,
               "Vitamins D and K given together")
})

# ---- Vitamin E (added 2026-11-04) -------------------------------------------
with_raw("Vitamin_E", "GSE110747", {   # mouse liver, vitamin E-enriched vs control diet x exercise, Mouse Gene 2.0 ST
  a <- array_dataset("Vitamin_E", "GSE110747", sym_fun = function(pl) parse_assignment(pl$gene_assignment))
  # deposited values are already log2 (e.g. Alb ~ 11.6), so they are used as is
  p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("VitE", p$title), "VitE_diet", "Control_diet"),
                   exercise = ifelse(grepl("Exercise", p$title), "Treadmill", "Sedentary"))
  save_dataset("GSE110747", "GSE110747", "Vitamin_E", "Mus musculus", "Affymetrix Mouse Gene 2.0 ST (GPL17400)",
               "log2 intensity (as deposited)", a$expr[, ph$sample], ph)
})

with_raw("Vitamin_E", "GSE1063", {   # rat, vitamin E-deficient vs sufficient diet, Clontech Atlas (~590 genes)
  a <- array_dataset("Vitamin_E", "GSE1063", sym_fun = function(pl) parse_paren_symbol(pl$GENE_NAME))
  p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("minus", p$title), "VitE_deficient", "VitE_sufficient"))
  save_dataset("GSE1063", "GSE1063", "Vitamin_E", "Rattus norvegicus", "Clontech Atlas Rat 1.2 (GPL155)",
               "log2 intensity", a$expr[, ph$sample], ph,
               "Small cDNA array; symbols parsed from gene names, so coverage is limited")
})

for (vE in list(c("GSE8150", "Affymetrix MOE430 2.0 (GPL1261)", "brain"), c("GSE8146", "Affymetrix MG-U74Av2 (GPL81)", "heart"))) {
  with_raw("Vitamin_E", vE[1], {   # old mice +/- tocopherol supplementation (and young controls), MAS5
    a <- array_dataset("Vitamin_E", vE[1]); p <- a$pheno
    arm <- ifelse(grepl("a- & g-tocopherol|agT", p$title), "Old_aT_gT",
           ifelse(grepl("a-tocopherol|aT-", p$title), "Old_aT",
           ifelse(grepl("young|C5", p$title), "Young_control", "Old_control")))
    ph <- data.frame(sample = p$gsm, condition = arm)
    save_dataset(vE[1], vE[1], "Vitamin_E", "Mus musculus", vE[2], "log2 MAS5", a$expr[, ph$sample], ph,
                 paste("Aged mouse", vE[3], "- tocopherol supplementation from middle age"))
  })
}

with_raw("Vitamin_E", "GSE195619", {   # HepaRG: untreated / oleic acid / oleic acid + three vitamin E sources
  a <- array_dataset("Vitamin_E", "GSE195619", sym_fun = function(pl) parse_assignment(pl$SPOT_ID))
  p <- a$pheno
  arm <- sub(".*\\[([A-Za-z-]+?)_?(OA)? ?[0-9]+ *\\]$", "\\1", p$title)
  cond <- ifelse(grepl("DMSO", p$title), "Untreated",
          ifelse(grepl("nVE", p$title), "OA_nVE", ifelse(grepl("sVE", p$title), "OA_sVE",
          ifelse(grepl("WGO", p$title), "OA_WGO_VE", "OA"))))
  ph <- data.frame(sample = p$gsm, condition = cond)
  save_dataset("GSE195619", "GSE195619", "Vitamin_E", "Homo sapiens", "Affymetrix GO_Screen_hu (GPL31262)",
               "log2 intensity", a$expr[, ph$sample], ph)
})

# ---- Second pass (added 2026-11-04): B1, B3, B5, iron, omega-3, vitamin C ------
sm_pheno <- function(group, acc) do.call(rbind, lapply(list.files(raw_path(group, acc, "matrix"), full.names = TRUE),
  function(f) { p <- read_series_matrix(f)$pheno; p[, c("gsm", "title", "source")] }))
ens_counts <- function(m, species = "Mus musculus") {   # Ensembl (versioned or not) rows -> symbol, log2 CPM
  ids <- sub("\\.\\d+$", "", rownames(m)); rownames(m) <- ids
  collapse_to_symbol(log_cpm(m), map_ids(ids, species, "ENSEMBL"))
}

with_raw("Vitamin_B1", "GSE31553", {   # human myotubes, benfotiamine vs none x normo/hyperglycaemia, 3 donors
  a <- array_dataset("Vitamin_B1", "GSE31553"); p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("benfotiamine", p$title), "Benfotiamine", "No_benfotiamine"),
                   glucose = ifelse(grepl("^hyper", p$title), "High", "Normal"), donor = sub(".*donor ", "", p$title))
  save_dataset("GSE31553", "GSE31553", "Vitamin_B1", "Homo sapiens", "Illumina HumanWG-6 v3 (GPL6884)",
               "log2 intensity", a$expr[, ph$sample], ph, "Benfotiamine (lipid-soluble thiamine) treatment, not deficiency")
})

with_raw("Vitamin_B3", "GSE345595", {   # ANDY mouse liver: niacin-free / niacin-deficient diets vs control (+ aging groups)
  x <- fread(raw_path("Vitamin_B3", "GSE345595", "suppl", "GSE345595_aging_liver_counts_matrix.csv.gz"), data.table = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x$Geneid; m <- ens_counts(m)
  p <- sm_pheno("Vitamin_B3", "GSE345595"); grp <- sub("^Liver_([^_]+)_.*", "\\1", p$title); smp <- sub(".*_", "", p$title)
  grp <- setNames(grp, smp)[colnames(m)]
  ph <- data.frame(sample = colnames(m), condition = unname(grp),
                   doxycycline = ifelse(grp %in% c("CD2W", "Old", "Young"), "No", "Yes"))
  save_dataset("GSE345595", "GSE345595", "Vitamin_B3", "Mus musculus", "RNA-seq", "log2 CPM", m, ph,
               "CD1D/CD2W = niacin-replete; ND1DL/ND1DS/ND2D = niacin-deficient; Old/Young = wild-type aging groups")
})

with_raw("Vitamin_B3", "GSE345344", {   # ANDY mouse liver: niacin-free diet, control, and recovery
  x <- fread(raw_path("Vitamin_B3", "GSE345344", "suppl", "GSE345344_recovery_liver_counts_matrix.csv.gz"), data.table = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x$Geneid; m <- ens_counts(m)
  p <- sm_pheno("Vitamin_B3", "GSE345344"); grp <- setNames(sub("^Liver_([^_]+)_.*", "\\1", p$title), sub(".*_", "", p$title))
  ph <- data.frame(sample = colnames(m), condition = unname(grp[colnames(m)]))
  save_dataset("GSE345344", "GSE345344", "Vitamin_B3", "Mus musculus", "RNA-seq", "log2 CPM", m, ph,
               "Same lab as GSE345595 but independent animals")
})

with_raw("Vitamin_B3", "GSE300498", {   # mouse intestine: control / methionine-restricted / tryptophan+niacin-free / both, 2 strains
  f <- raw_path("Vitamin_B3", "GSE300498", "suppl", "GSE300498_processed_data.txt.gz")
  x <- fread(f, skip = 1, header = TRUE, data.table = FALSE)
  x <- x[, c(1, setdiff(which(!grepl("FMT", names(x))), 1))]   # faecal-transfer columns are mostly NA / near-empty
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  colnames(m) <- make.unique(gsub("[ /]", "_", sub(" mice$", "", colnames(m))))
  m <- ens_counts(m)
  arm <- sub("_.*", "", colnames(m))
  ph <- data.frame(sample = colnames(m), condition = arm, strain = sub("^[^_]+_", "", colnames(m)))
  save_dataset("GSE300498", "GSE300498", "Vitamin_B3", "Mus musculus", "RNA-seq", "log2 CPM", m, ph,
               "TNR/MTNR lack tryptophan and niacin; r* = 2-day recovery; faecal-transfer samples omitted (empty columns); one sample per strain per arm")
})

with_raw("Vitamin_B3", "GSE92499", {   # mouse liver: nicotinamide none / low / high x standard / high-fat diet (Illumina)
  x <- fread(raw_path("Vitamin_B3", "GSE92499", "suppl", "GSE92499_non-normalized.txt.gz"), skip = "ID\t", data.table = FALSE)
  keep <- names(x) != "ID" & !grepl("Detection|Pval", names(x)) & vapply(x, is.numeric, TRUE)
  m <- as.matrix(x[, keep]); rownames(m) <- x$ID
  m <- limma::normalizeBetweenArrays(log2(pmax(m, 1)), method = "quantile")
  pl <- read_soft_platform(raw_path("Vitamin_B3", "GSE92499", "soft", "GSE92499_family.soft.gz"))
  m <- collapse_to_symbol(m, pl[[symbol_col(pl)]][match(rownames(m), pl$ID)])
  s <- colnames(m)
  ph <- data.frame(sample = s,
                   condition = ifelse(grepl("-H-", s, ignore.case = TRUE), "High_NAM", ifelse(grepl("-L-", s, ignore.case = TRUE), "Low_NAM", "No_NAM")),
                   diet = ifelse(startsWith(s, "HFD"), "HFD", "SD"))
  save_dataset("GSE92499", "GSE92499", "Vitamin_B3", "Mus musculus", "Illumina MouseWG-6 v2 (GPL6885)",
               "log2 intensity (quantile-normalised from non-normalised file)", m, ph)
})

with_raw("Vitamin_B5", "GSE130705", {   # mouse brain, Pank1/Pank2 double knockout vs WT, day 10/20, both sexes (STAR counts)
  tarf <- raw_path("Vitamin_B5", "GSE130705", "suppl", "GSE130705_RAW.tar"); td <- tempfile(); dir.create(td); untar(tarf, exdir = td)
  p <- read_series_matrix(list.files(raw_path("Vitamin_B5", "GSE130705", "matrix"), full.names = TRUE)[1])$pheno
  per <- lapply(p$gsm, function(gsm) { y <- fread(list.files(td, pattern = paste0("^", gsm, "_"), full.names = TRUE), header = FALSE, data.table = FALSE)
    y <- y[!startsWith(y[[1]], "N_"), ]; setNames(y[[2]], y[[1]]) })
  unlink(td, recursive = TRUE); g <- Reduce(intersect, lapply(per, names))
  m <- vapply(per, function(v) v[g], numeric(length(g))); rownames(m) <- g; colnames(m) <- p$gsm
  m <- ens_counts(m)
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("knockout", p$genotype.variation), "Pank1_Pank2_dKO", "WT"),
                   sex = p$Sex, day = paste0("day", p$day))
  save_dataset("GSE130705", "GSE130705", "Vitamin_B5", "Mus musculus", "RNA-seq (GPL17021)", "log2 CPM", m, ph,
               "Genetic coenzyme A deficiency (pantothenate kinase knockout), not dietary")
})

with_raw("Vitamin_B5", "GSE187456", {   # mouse Tc0 CD8 T cells +/- coenzyme A (Tc22 cells excluded)
  x <- fread(raw_path("Vitamin_B5", "GSE187456", "suppl", "GSE187456_stpaul_saibil_ohashi_2021_counts.txt.gz"), data.table = FALSE)
  m <- as.matrix(x[, grep("^tc0_", names(x))]); rownames(m) <- x$ensembl_gene_id; m <- ens_counts(m)
  ph <- data.frame(sample = colnames(m), condition = ifelse(grepl("_coa_", colnames(m)), "CoA", "Untreated"),
                   mouse = sub(".*_", "", colnames(m)))
  save_dataset("GSE187456", "GSE187456", "Vitamin_B5", "Mus musculus", "RNA-seq (GPL13112)", "log2 CPM", m, ph)
})

with_raw("Iron", "GSE162493", {   # mouse heart, iron-deficient (<2 ppm) vs control diet (featureCounts, gene symbols)
  x <- fread(raw_path("Iron", "GSE162493", "suppl", "GSE162493_featuresCount.txt.gz"), skip = "Geneid", data.table = FALSE)
  m <- as.matrix(x[, grep("bam$", names(x))]); colnames(m) <- sub("_.*", "", colnames(m))
  m <- collapse_to_symbol(log_cpm(m), x$Geneid)
  ph <- data.frame(sample = colnames(m), condition = ifelse(startsWith(colnames(m), "ID"), "Iron_deficient", "Control"))
  save_dataset("GSE162493", "GSE162493", "Iron", "Mus musculus", "RNA-seq (GPL21103)", "log2 CPM", m, ph)
})

with_raw("Iron", "GSE157667", {   # mouse placenta E12.5, iron-deficient vs adequate, both sexes
  x <- fread(raw_path("Iron", "GSE157667", "suppl", "GSE157667_featurecounts_mus_rm.txt.gz"), data.table = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x$Geneid; m <- ens_counts(m)
  s <- colnames(m)   # GEO descriptions: D = female ID, MD = male ID, A = female IA, MA = male IA
  ph <- data.frame(sample = s, condition = ifelse(grepl("D", sub("\\d+$", "", s)), "Iron_deficient", "Iron_adequate"),
                   sex = ifelse(startsWith(s, "M"), "Male", "Female"))
  save_dataset("GSE157667", "GSE157667", "Iron", "Mus musculus", "RNA-seq (GPL21103)", "log2 CPM", m, ph)
})

with_raw("Iron", "GSE30533", {   # rat liver, iron-deficient vs control diet, Affymetrix Rat 230 2.0
  a <- array_dataset("Iron", "GSE30533"); p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("iron-deficient", p$title), "Iron_deficient", "Control"))
  save_dataset("GSE30533", "GSE30533", "Iron", "Rattus norvegicus", "Affymetrix Rat 230 2.0 (GPL1355)", "log2 intensity",
               a$expr[, ph$sample], ph)
})

with_raw("Iron", "GSE175739", {   # heart-failure patients with vs without systemic iron deficiency, whole blood, matched pairs
  # The GEO platform table gives accessions for only ~1/3 of transcript clusters (TFRC is missing),
  # so symbols come from the Bioconductor annotation package for this array
  if (!requireNamespace("hugene20sttranscriptcluster.db", quietly = TRUE))
    stop("BiocManager::install('hugene20sttranscriptcluster.db') is needed for GSE175739")
  hg_symbol <- function(pl) suppressMessages(AnnotationDbi::mapIds(
    hugene20sttranscriptcluster.db::hugene20sttranscriptcluster.db, as.character(pl$ID), "SYMBOL", "PROBEID", multiVals = "first"))
  a <- array_dataset("Iron", "GSE175739", sym_fun = hg_symbol); p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(p$iron.status == "1", "Iron_deficient", "Iron_replete"),
                   pair = paste0("pair", p$sample.pair.id))
  save_dataset("GSE175739", "GSE175739", "Iron", "Homo sapiens", "Affymetrix HuGene 2.0 ST (GPL16686)", "log2 intensity",
               a$expr[, ph$sample], ph, "Patients (heart failure) - pairs matched on sex, age, Hb, LVEF, NT-proBNP, aetiology, CrCl")
})

with_raw("Iron", "GSE314015", {   # human arterial and lymphatic endothelial cells, deferoxamine vs vehicle
  x <- fread(raw_path("Iron", "GSE314015", "suppl", "GSE314015_counts.csv.gz"), data.table = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m <- ens_counts(m, "Homo sapiens")
  s <- colnames(m)
  ph <- data.frame(sample = s, condition = ifelse(substr(s, 2, 2) == "D", "Deferoxamine", "Vehicle"),
                   cell_type = ifelse(startsWith(s, "A"), "Coronary_artery_EC", "Lymphatic_EC"))
  save_dataset("GSE314015", "GSE314015", "Iron", "Homo sapiens", "RNA-seq", "log2 CPM", m, ph,
               "Iron chelation in culture, not dietary deficiency")
})

with_raw("Omega_3", "GSE48368", {   # human PBMC, fish oil vs high-oleic sunflower oil, weeks 0/3/7
  a <- array_dataset("Omega_3", "GSE48368"); p <- a$pheno
  oil <- ifelse(grepl("fish", p$treatment), "FO", "HOSO"); wk <- sub("wk ", "wk", p$time.point)
  ph <- data.frame(sample = p$gsm, condition = paste(oil, wk, sep = "_"), week = wk, subject = p$individual.id)
  save_dataset("GSE48368", "GSE48368", "Omega_3", "Homo sapiens", "Illumina HumanHT-12 v4 (GPL10558)", "log2 intensity",
               a$expr[, ph$sample], ph)
})

with_raw("Omega_3", "GSE29680", {   # Alzheimer patients' PBMC, DHA-rich n-3 vs placebo, before/after 6 months
  a <- array_dataset("Omega_3", "GSE29680"); p <- a$pheno
  arm <- ifelse(grepl("placebo", p$title), "Placebo", "n3"); t <- ifelse(grepl("before", p$title), "before", "after")
  ph <- data.frame(sample = p$gsm, condition = paste(arm, t, sep = "_"))
  save_dataset("GSE29680", "GSE29680", "Omega_3", "Homo sapiens", "Affymetrix HG-Focus (GPL201)", "log2 MAS5",
               a$expr[, ph$sample], ph, "Patients (Alzheimer's disease)")
})

with_raw("Vitamin_C", "GSE295923", {   # E8.5 headfolds of Gulo-/- embryos, maternal vitamin C minus vs plus
  tarf <- raw_path("Vitamin_C", "GSE295923", "suppl", "GSE295923_RAW.tar"); td <- tempfile(); dir.create(td)
  fs <- grep("countMatrix_by_featureCounts", untar(tarf, list = TRUE), value = TRUE); untar(tarf, files = fs, exdir = td)
  per <- lapply(fs, function(f) { y <- fread(file.path(td, f), skip = "Geneid", data.table = FALSE); setNames(y[[8]], y$gene_name) })
  unlink(td, recursive = TRUE)
  g <- Reduce(intersect, lapply(per, names))
  m <- vapply(per, function(v) tapply(v, names(v), sum)[g], numeric(length(g)))
  nm <- sub("\\.countMatrix.*", "", sub("^GSM\\d+_RNA_", "", basename(fs))); colnames(m) <- nm; rownames(m) <- g
  m <- collapse_to_symbol(log_cpm(m), g)
  ph <- data.frame(sample = colnames(m), condition = ifelse(grepl("AAmin", colnames(m)), "VitC_minus", "VitC_plus"),
                   background = sub("_.*", "", colnames(m)), sex = ifelse(grepl("female", colnames(m)), "Female", "Male"))
  save_dataset("GSE295923", "GSE295923", "Vitamin_C", "Mus musculus", "RNA-seq (GPL24247)", "log2 CPM", m, ph,
               "Maternal vitamin C status of Gulo-/- dams; embryonic headfold")
})

# ---- Index ------------------------------------------------------------------
idx <- do.call(rbind, index)
old_idx_file <- file.path(OUT, "dataset_index.csv")
if (file.exists(old_idx_file)) {               # keep rows for datasets not re-processed this run
  old_idx <- fread(old_idx_file, data.table = FALSE)
  keep <- old_idx[old_idx$accession %in% skipped_acc & !old_idx$dataset %in% idx$dataset, names(idx)]
  idx <- rbind(keep, idx)
}
fwrite(idx, file.path(OUT, "dataset_index.csv"))
message("==== Done: ", nrow(idx), " datasets. Index: ", file.path(OUT, "dataset_index.csv"))
