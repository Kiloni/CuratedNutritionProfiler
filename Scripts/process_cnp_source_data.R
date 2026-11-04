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

array_dataset <- function(group, acc, symbol_from = c("platform", "entrez")) {
  symbol_from <- match.arg(symbol_from)
  sm <- read_series_matrix(list.files(raw_path(group, acc, "matrix"), full.names = TRUE)[1])
  x  <- log_array(sm$expr)
  sym <- if (symbol_from == "platform") {
    pl <- read_soft_platform(raw_path(group, acc, "soft", paste0(acc, "_family.soft.gz")))
    pl[[symbol_col(pl)]][match(rownames(x), as.character(pl$ID))]
  } else {
    map_ids(sub("_at$", "", rownames(x)), "Homo sapiens", "ENTREZID")
  }
  list(expr = collapse_to_symbol(x, sym), pheno = sm$pheno)
}

ncbi_counts_dataset <- function(group, acc) {
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
local({
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
local({
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
local({
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
local({
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
local({   # GSE63086: MG63 cells, vehicle vs 1,25D x control vs HNRNPC knock-down
  d <- ncbi_counts_dataset("Vitamin_D", "GSE63086"); p <- d$pheno
  ph <- data.frame(sample = p$gsm,
                   condition = ifelse(grepl("^Vitamin D", p$title), "Calcitriol", "Vehicle"),
                   knockdown = ifelse(grepl("KD", p$title), "HNRNPC_KD", "Control"))
  save_dataset("GSE63086", "GSE63086", "Vitamin_D", "Homo sapiens", "RNA-seq (GPL16791)",
               "log2 CPM (NCBI counts)", d$expr[, intersect(ph$sample, colnames(d$expr))], ph)
})

local({   # GSE221690: 4 prostate cell lines, vehicle vs 1,25D3
  d <- ncbi_counts_dataset("Vitamin_D", "GSE221690"); p <- d$pheno
  ph <- data.frame(sample = p$gsm,
                   condition = ifelse(grepl("treat", p$title), "Calcitriol", "Vehicle"),
                   cell_line = p$cell.line)
  save_dataset("GSE221690", "GSE221690", "Vitamin_D", "Homo sapiens", "RNA-seq (GPL16791)",
               "log2 CPM (NCBI counts)", d$expr[, intersect(ph$sample, colnames(d$expr))], ph)
})

local({   # GSE15947: RWPE1 cells, vehicle vs 1,25D at 6/24/48 h, GPL570 (already log2)
  a <- array_dataset("Vitamin_D", "GSE15947"); p <- a$pheno
  ph <- data.frame(sample = p$gsm,
                   condition = ifelse(grepl("Vehicle", p$treatment), "Vehicle", "Calcitriol"),
                   time = p$time.point)
  save_dataset("GSE15947", "GSE15947", "Vitamin_D", "Homo sapiens",
               "Affymetrix HG-U133 Plus 2 (GPL570)", "log2 intensity (RMA)", a$expr[, ph$sample], ph)
})

local({   # GSE283020: mouse TH2 cells, regular vs vitamin D-deficient diet, counts keyed by Ensembl
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
local({   # GSE233598: mouse liver, Gulo-/- vitamin C dose series + WT, TPM
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

local({   # GSE190371: mouse kidney scRNA-seq, Gulo-/- control vs vitamin C-deficient -> pseudobulk
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
local({   # GSE161763: mouse embryonic brain, Ronin F80L (cblX) vs control, counts
  x <- fread(raw_path("Vitamin_B12", "GSE161763", "suppl", "GSE161763_RONIN_F80L_RNAseq_raw_counts.txt.gz"),
             data.table = FALSE)
  g <- x$Gene; x <- as.matrix(x[, -1]); rownames(x) <- g
  m <- collapse_to_symbol(log_cpm(x), g)
  ph <- data.frame(sample = colnames(m), condition = ifelse(grepl("mutant", colnames(m)), "cblX_F80L", "Control"))
  save_dataset("GSE161763", "GSE161763", "Vitamin_B12", "Mus musculus", "RNA-seq (GPL17021)", "log2 CPM", m, ph)
})

local({   # GSE103417: N1E-115 cells, WT / OT (no Cbl binding) / TO (sequesters Cbl), Agilent
  a <- array_dataset("Vitamin_B12", "GSE103417"); p <- a$pheno
  arm <- sub("^(proliferated|differentiated)_neuroblastoma_", "", p$source)
  ph <- data.frame(sample = p$gsm,
                   condition = c(WT = "WT", OT = "OT_control_fusion", TO = "TO_B12_sequestered")[arm],
                   medium = p$medium)
  save_dataset("GSE103417", "GSE103417", "Vitamin_B12", "Mus musculus",
               "Agilent Whole Mouse Genome (GPL4134)", "log2 intensity", a$expr[, ph$sample], ph)
})

local({   # GSE104164: rat cerebellum, control vs methyl-donor (B12/folate) deficiency, Agilent
  a <- array_dataset("Vitamin_B12", "GSE104164"); p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = ifelse(grepl("^MDD", p$source), "MDD", "Control"))
  save_dataset("GSE104164", "GSE104164", "Vitamin_B12", "Rattus norvegicus",
               "Agilent Rat Gene Expression (GPL14797)", "log2 intensity", a$expr[, ph$sample], ph)
})

local({   # GSE99113: mouse astrocytes, B12-sufficient vs B12-free medium (2 vs 2)
  x <- fread(raw_path("Vitamin_B12", "GSE99113", "suppl", "GSE99113_JC_NS190_AstrocyteB12.txt.gz"),
             data.table = FALSE)
  g <- x$name; x <- as.matrix(x[, -1]); rownames(x) <- g
  m <- collapse_to_symbol(log_plus(x), g)
  ph <- data.frame(sample = colnames(m), condition = ifelse(grepl("def", colnames(m)), "B12_deficient", "B12_sufficient"))
  save_dataset("GSE99113", "GSE99113", "Vitamin_B12", "Mus musculus", "RNA-seq (GPL19057)",
               "log2(normalised expression+1)", m, ph, "2 vs 2 - descriptive only")
})

# ---- Vitamin K --------------------------------------------------------------
local({   # GSE148011: p53/PTEN-null MEFs, vitamin K precursor (MSB) 0/10/20 uM (+ erastin arms)
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
local({   # GSE12375: human PBMC, EPA/DHA vs MUFA (oleic acid) at baseline and 6 months, NuGO array
  a <- array_dataset("Omega_3", "GSE12375", symbol_from = "entrez"); p <- a$pheno
  ph <- data.frame(sample = p$gsm, condition = sub(".*_(EPADHA|MUFA)_(BASE|END)$", "\\1_\\2", p$title),
                   subject = sub("^A101_([^_]+)_.*", "\\1", p$title))
  save_dataset("GSE12375", "GSE12375", "Omega_3", "Homo sapiens",
               "NuGO Hs1a520180 custom Affymetrix (GPL7144)", "log2 intensity", a$expr[, ph$sample], ph)
})

# ---- Index ------------------------------------------------------------------
idx <- do.call(rbind, index)
fwrite(idx, file.path(OUT, "dataset_index.csv"))
message("==== Done: ", nrow(idx), " datasets. Index: ", file.path(OUT, "dataset_index.csv"))
