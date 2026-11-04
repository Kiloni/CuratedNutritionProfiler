# =============================================================================
# CuratedNutritionProfiler (CNP) — download fresh copies of all source data
# behind the nutrition signatures, for re-analysis.
#
# What it does
#   1. GEO: for every series (GSE) below, downloads the series matrix file(s),
#      the SOFT family file, every supplementary file, and — when GEO has
#      generated them — NCBI's uniformly processed RNA-seq counts (raw counts,
#      TPM, FPKM + gene annotation). GDS4860 is fetched as a GDS SOFT file.
#   2. GWAS Catalog (REST API v2): all associations for each study (GCST) and
#      each trait (EFO/OBA) used to build the GWAS-based signatures, saved as
#      raw JSON plus a flat TSV, together with study metadata.
#   3. Writes a manifest (accession, file, URL, size, md5, status) and the R
#      session info, so the downloads can be cited in the methods section.
#
# Safe to re-run: files already downloaded at the same size are skipped.
# Nothing outside OUT_DIR is touched.
#
# Needs: R >= 4.1 and the 'jsonlite' package (installed below if missing).
# Time/space: ~3.3 GB in total (survey 2026-11-04); files above MAX_FILE_GB
#   are skipped and listed in the manifest so you can decide on them
#   separately. As of 2026-11-04 the only such file is GSE271249_RAW.tar
#   (6.6 GB, Retinoic acid:YAP1 candidate) — deferred until needed.
#
# Notes from the 2026-11-04 review:
#   - Output now lives inside this R project (data-raw/CNP_Source_Data/<group>/)
#     until the nutrition signatures list is validated; the folder is
#     git-ignored.
#   - EFO_0004857 (retinol measurement) has been retired from the GWAS
#     Catalog; its associations are now under OBA_1000995 ("amount of
#     retinol"). Both are queried so the manifest records the change.
#   - Many GCST9xxxxxxx studies are summary-statistics-only and return 0
#     curated associations. That is expected; the study metadata JSON holds
#     the full_summary_stats FTP link if those are needed later.
#   - The manifest is now written after every accession, and GWAS requests
#     are wrapped in tryCatch, so a failure part-way through no longer loses
#     the record of what was already downloaded.
# =============================================================================

# ---- Settings ---------------------------------------------------------------
PROJECT_DIR <- if (file.exists("curatedNutritionDataset.Rproj")) getwd() else
  path.expand("~/RStudio Projects/curatedNutritionDataset")
OUT_DIR     <- file.path(PROJECT_DIR, "data-raw", "CNP_Source_Data")
MAX_FILE_GB <- 5          # skip any single file larger than this (set Inf to take all)
DOWNLOAD_NCBI_COUNTS <- TRUE
# Skip GEO accessions already converted by Scripts/process_cnp_source_data.R (their raw
# files are deleted after processing to save space). Set FALSE to re-download everything.
SKIP_PROCESSED <- TRUE
# Accessions reviewed on 2026-11-04 and excluded from signature validation (no usable nutrient
# contrast, or not expression data); their raw files were deleted. Remove an ID here to fetch it.
#   GSE266352 LAMN vs CRC metastases | GSE124395 liver cell atlas | GSE246328 NASH diet, no iron arm
#   GSE296260 ccRCC vs kidney | GSE143307 DNA methylation (vitamin K responders)
#   GSE271214 / GSE271249 RA-treated hESC, 2 per arm (RAW tar 7.1 GB) | GDS4860 subset of GSE46268
#   GSE19292 CodeLink rat array without gene symbols, 2 vs 2
EXCLUDED <- c("GSE266352", "GSE124395", "GSE246328", "GSE296260", "GSE143307", "GSE271214", "GSE271249", "GDS4860",
              "GSE19292")

options(timeout = 3600)   # allow up to 1 h per file
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(jsonlite)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
manifest_file <- file.path(OUT_DIR, "download_manifest.csv")
log_file      <- file.path(OUT_DIR, "download_log.txt")
processed_index <- file.path(PROJECT_DIR, "data", "processed", "dataset_index.csv")
processed_acc <- if (SKIP_PROCESSED && file.exists(processed_index))
  unique(read.csv(processed_index, stringsAsFactors = FALSE)$accession) else character(0)

# ---- Accessions -------------------------------------------------------------
# group = folder the files go in; role = why the dataset is here.
geo <- read.csv(text = '
group,accession,role
Vitamin_A,GSE134933,"A_Retinol_GLM_43 training (VAD vs VAS mouse kidney)"
Vitamin_A,GSE299530,"A_Retinol_Quiles_20 (VAD mouse lung)"
Vitamin_A,GSE266352,"Retinol human training (Retinol_Training_SE.R)"
Vitamin_A,GDS4860,"Retinol scoring dataset (GEO DataSet)"
Vitamin_A,GSE46268,"Parent series of GDS4860"
Vitamin_A,GSE37313,"Retinol_GSE37313 scoring dataset"
Vitamin_A,GSE271249,"Retinoic acid:YAP1 (metadata in Gene Signatures) - candidate"
Vitamin_A,GSE271214,"Retinoic acid:YAP1 scRNA-seq - candidate"
Vitamin_B12,GSE103417,"B12.R input"
Vitamin_B12,GSE104164,"B12.R input"
Vitamin_B12,GSE161763,"B12.R input"
Vitamin_B12,GSE99113,"B12.R input"
Vitamin_C,GSE233598,"Vitamin_C signature source"
Vitamin_C,GSE190371,"Vitamin C deficiency - candidate"
Vitamin_D,GSE63086,"Vitamin D human training (VitaminD_Training_SE.R)"
Vitamin_D,GSE283020,"Vitamin D mouse TH2 training (VitaminD_Training_SE.R)"
Vitamin_D,GSE15947,"1,25(OH)2D RWPE1 time course - candidate"
Vitamin_D,GSE221690,"D3 vs EtOH RC43 cells - candidate"
Vitamin_D,GSE124395,"Human liver cell atlas (Aizarani MSigDB sets) - candidate"
Vitamin_K,GSE148011,"Vitamin_K signature source"
Vitamin_K,GSE143307,"Accession list entry for K (listed as GSA143307)"
Iron,GSE246328,"Iron_Homeostasis testing"
Omega_3,GSE12375,"Omega-3 fish-oil PBMC testing"
Unverified,GSE296260,"In Downloads, nutrient not confirmed"
Vitamin_B6,GSE345689,"Low vs high vitamin B6 diet, mouse kidney (WT and Sgpl1 R222Q) - added 2026-11-04"
Vitamin_B6,GSE139892,"PDXK inhibitor 4-O-methylpyridoxine vs vehicle, mouse leukemic cells - added 2026-11-04"
Folate,GSE275393,"Maternal folic acid 0.1 / 3 / 30 ppm x Tet1 genotype, E11.5 brain - added 2026-11-04"
Folate,GSE324232,"Folate-deficient vs sufficient spermatocytes/spermatids - added 2026-11-04"
Folate,GSE271941,"Moderate vs fortified folic acid MEFs x Rsg1 genotype - added 2026-11-04"
Carnitine,GSE343742,"L-carnitine vs none, +/- LPS, mouse BMM - added 2026-11-04"
Carnitine,GSE344237,"Resveratrol + L-carnitine vs none x diet, mouse liver - added 2026-11-04"
Vitamin_K,GSE267133,"Erastin +/- vitamin K1, ARPE-19 cells - added 2026-11-04"
Vitamin_K,GSE237219,"Vitamin K1 vs control, Jurkat cells - added 2026-11-04"
Vitamin_K,GSE272401,"RANKL +/- vitamins D and K, mouse BMDM - added 2026-11-04"
Vitamin_E,GSE1063,"Vitamin E deficient vs sufficient diet, rat - added 2026-11-04"
Vitamin_E,GSE110747,"Vitamin E-enriched vs control diet x exercise, mouse liver - added 2026-11-04"
Vitamin_E,GSE8150,"Old mice +/- tocopherol supplementation, brain - added 2026-11-04"
Vitamin_E,GSE8146,"Old mice +/- tocopherol supplementation, heart - added 2026-11-04"
Vitamin_E,GSE195619,"HepaRG oleic acid +/- vitamin E sources - added 2026-11-04"
Vitamin_B1,GSE31553,"Benfotiamine vs none, human myotubes x glucose - added 2026-11-04 (2nd pass)"
Vitamin_B1,GSE19292,"Thiamine supplementation vs none, OLETF rat liver - added 2026-11-04 (2nd pass)"
Vitamin_B3,GSE345595,"Niacin-free vs control diet, ANDY mouse liver - added 2026-11-04 (2nd pass)"
Vitamin_B3,GSE345344,"Niacin-free diet and recovery, ANDY mouse liver - added 2026-11-04 (2nd pass)"
Vitamin_B3,GSE300498,"Tryptophan/niacin-free and methionine-restricted diets, mouse intestine - added 2026-11-04 (2nd pass)"
Vitamin_B3,GSE92499,"Nicotinamide supplementation x diet, mouse liver - added 2026-11-04 (2nd pass)"
Vitamin_B5,GSE130705,"Pank1/Pank2 double knockout (CoA deficiency), mouse brain - added 2026-11-04 (2nd pass)"
Vitamin_B5,GSE187456,"Coenzyme A vs untreated, mouse Tc0 CD8 T cells - added 2026-11-04 (2nd pass)"
Iron,GSE162493,"Iron-deficient vs control diet, mouse heart - added 2026-11-04 (2nd pass)"
Iron,GSE157667,"Iron-deficient vs adequate, mouse placenta - added 2026-11-04 (2nd pass)"
Iron,GSE30533,"Iron-deficient vs control diet, rat liver - added 2026-11-04 (2nd pass)"
Iron,GSE175739,"Heart failure patients with vs without iron deficiency, whole blood - added 2026-11-04 (2nd pass)"
Iron,GSE314015,"Deferoxamine vs vehicle, human endothelial cells - added 2026-11-04 (2nd pass)"
Omega_3,GSE48368,"Fish oil vs high-oleic sunflower oil, human PBMC wk 0/3/7 - added 2026-11-04 (2nd pass)"
Omega_3,GSE29680,"DHA-rich n-3 vs placebo, Alzheimer patient PBMC - added 2026-11-04 (2nd pass)"
Vitamin_C,GSE295923,"Gulo-/- embryos +/- maternal vitamin C, E8.5 headfold - added 2026-11-04 (2nd pass)"
', stringsAsFactors = FALSE)

gwas_studies <- read.csv(text = '
group,accession,role
Vitamin_A,GCST90265907,"Retinol GWAS (accession list)"
Vitamin_B1,GCST90300747,"Thiamine GWAS (accession list)"
Vitamin_B2,GCST90255448,"Riboflavin/Niacin GWAS (accession list)"
Vitamin_B3,GCST90300148,"Niacinamide GWAS (accession list)"
Vitamin_B5,GCST90013572,"Pantothenic acid GWAS (accession list)"
Vitamin_B5,GCST90176068,"Pantothenic acid GWAS (accession list)"
Vitamin_B6,GCST90257552,"Pyridoxine GWAS (accession list)"
Vitamin_B6,GCST90258157,"Pyridoxine GWAS (accession list)"
Vitamin_B6_B9_B12,GCST002559,"B6/B9/B12 GWAS (accession list)"
Vitamin_B12,GCST90448636,"B12 GWAS (accession list)"
Vitamin_B12,GCST90277442,"B12 GWAS (accession list)"
Vitamin_B12,GCST001762,"B12 GWAS (accession list)"
GWAS_other,GCST001216,"Listed in Nutritional_Signatures.numbers"
GWAS_other,GCST009391,"Listed in Nutritional_Signatures.numbers"
GWAS_other,GCST90139670,"Listed in Nutritional_Signatures.numbers"
', stringsAsFactors = FALSE)

gwas_traits <- read.csv(text = '
group,accession,role
Vitamin_A,EFO_0004857,"Retinol_GWAS_EFO_0004857 export used for A_Retinol_GWAS (retired ID; kept for the record)"
Vitamin_A,OBA_1000995,"Current GWAS Catalog ID for retinol (amount of retinol); replaces EFO_0004857"
Vitamin_B1,EFO_0010540,"Thiamine_GWAS_EFO_0010540 export used for B1_Thiamine"
Vitamin_B3,OBA_2040014,"Vitamin-B3 GWAS export used for B3_Niacin"
Vitamin_B5,EFO_0010519,"Vitamin-B5 GWAS export used for B5_Pantothenic Acid"
Vitamin_B6,EFO_0004621,"Vitamin B6 measurement - source of the B6_Pyridoxine gene list"
Vitamin_E,EFO_0007898,"Alpha-tocopherol measurement - basis for Vitamin_E_GWAS"
Vitamin_E,EFO_0004738,"Tocopherol measurement - basis for Vitamin_E_GWAS"
Vitamin_K,EFO_0004618,"Vitamin K measurement - basis for Vitamin_K_GWAS"
Carnitine,EFO_0010469,"Carnitine measurement - basis for Carnitine_GWAS"
', stringsAsFactors = FALSE)

# ---- Helpers ----------------------------------------------------------------
logmsg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(...))
  message(msg); cat(msg, "\n", file = log_file, append = TRUE)
}

manifest <- data.frame()
add_manifest <- function(group, accession, file, url, status, role = "") {
  size <- if (file.exists(file)) file.size(file) else NA
  md5  <- if (file.exists(file)) unname(tools::md5sum(file)) else NA
  manifest <<- rbind(manifest, data.frame(
    group = group, accession = accession, role = role,
    file = sub(paste0(OUT_DIR, "/"), "", file, fixed = TRUE), url = url,
    bytes = size, md5 = md5, status = status,
    downloaded = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE))
}
# Rows from earlier runs are kept for every accession this run does not touch
old_manifest <- if (file.exists(manifest_file)) read.csv(manifest_file, stringsAsFactors = FALSE) else NULL
save_manifest <- function() {
  if (!nrow(manifest)) return(invisible())
  m <- manifest
  if (!is.null(old_manifest)) m <- rbind(old_manifest[!old_manifest$accession %in% manifest$accession, names(manifest)], manifest)
  write.csv(m, manifest_file, row.names = FALSE)
}

remote_size <- function(url) {
  h <- tryCatch(curlGetHeaders(url, redirect = TRUE), error = function(e) NULL)
  if (is.null(h)) return(list(ok = FALSE, bytes = NA))
  st <- attr(h, "status")
  cl <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
  bytes <- if (length(cl)) as.numeric(trimws(sub("^[^:]+:", "", tail(cl, 1)))) else NA
  list(ok = !is.null(st) && st < 400, bytes = bytes)
}

fetch <- function(url, dest, group, accession, role = "", quiet_404 = FALSE) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  rs <- remote_size(url)
  if (!rs$ok) {
    if (!quiet_404) logmsg("  not available: ", url)
    if (!quiet_404) add_manifest(group, accession, dest, url, "not available", role)
    return(invisible(FALSE))
  }
  if (!is.na(rs$bytes) && rs$bytes > MAX_FILE_GB * 1e9) {
    logmsg("  SKIPPED (", round(rs$bytes / 1e9, 2), " GB > cap): ", url)
    add_manifest(group, accession, dest, url,
                 sprintf("skipped: %.2f GB over cap", rs$bytes / 1e9), role)
    return(invisible(FALSE))
  }
  if (file.exists(dest) && !is.na(rs$bytes) && file.size(dest) == rs$bytes) {
    add_manifest(group, accession, dest, url, "already present", role)
    return(invisible(TRUE))
  }
  ok <- tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl"); TRUE
  }, error = function(e) { logmsg("  FAILED: ", url, " (", conditionMessage(e), ")"); FALSE })
  # A .gz that is not really gzip is an error page from the server, not data
  if (ok && grepl("\\.gz$", dest)) {
    magic <- readBin(dest, "raw", 2)
    if (!identical(magic, as.raw(c(0x1f, 0x8b)))) {
      unlink(dest); ok <- FALSE
      if (!quiet_404) logmsg("  not a valid gzip (server error page): ", url)
      if (quiet_404) return(invisible(FALSE))
    }
  }
  if (ok) logmsg("  ok  ", basename(dest), "  ", format(structure(file.size(dest), class = "object_size"), units = "auto"))
  add_manifest(group, accession, dest, url, if (ok) "downloaded" else "failed", role)
  invisible(ok)
}

# Files listed in a GEO FTP directory (empty if the directory does not exist)
list_ftp_dir <- function(url) {
  html <- tryCatch(readLines(url, warn = FALSE), error = function(e) character(0))
  hrefs <- regmatches(html, regexpr('href="[^"]+"', html))
  hrefs <- gsub('^href="|"$', "", hrefs)
  hrefs[!grepl("/$|^\\?|^/|^https?:", hrefs)]
}

geo_stub <- function(acc) {        # GSE103417 -> GSE103nnn ; GSE7451 -> GSE7nnn ; GSE123 -> GSEnnn
  prefix <- sub("\\d+$", "", acc); num <- sub("^[A-Z]+", "", acc)
  paste0(prefix, if (nchar(num) <= 3) "" else substr(num, 1, nchar(num) - 3), "nnn")
}

download_gse <- function(acc, group, role) {
  logmsg(acc, " (", group, ")")
  base <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s", geo_stub(acc), acc)
  dest_dir <- file.path(OUT_DIR, group, acc)

  # Series matrix file(s) — one per platform
  for (f in list_ftp_dir(paste0(base, "/matrix/")))
    fetch(paste0(base, "/matrix/", f), file.path(dest_dir, "matrix", f), group, acc, role)

  # SOFT family file (sample-level metadata for every GSM)
  fetch(sprintf("%s/soft/%s_family.soft.gz", base, acc),
        file.path(dest_dir, "soft", paste0(acc, "_family.soft.gz")), group, acc, role)

  # All supplementary files deposited by the authors
  supp <- list_ftp_dir(paste0(base, "/suppl/"))
  if (!length(supp)) logmsg("  no supplementary files listed")
  for (f in supp)
    fetch(paste0(base, "/suppl/", f), file.path(dest_dir, "suppl", f), group, acc, role)

  # NCBI-generated RNA-seq counts (only exist for some human/mouse RNA-seq series)
  if (DOWNLOAD_NCBI_COUNTS) {
    gen <- c(
      sprintf("%s_raw_counts_GRCh38.p13_NCBI.tsv.gz", acc),
      sprintf("%s_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz", acc),
      sprintf("%s_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz", acc),
      sprintf("%s_raw_counts_GRCm39_NCBI.tsv.gz", acc),
      sprintf("%s_norm_counts_TPM_GRCm39_NCBI.tsv.gz", acc),
      sprintf("%s_norm_counts_FPKM_GRCm39_NCBI.tsv.gz", acc))
    # NCBI's gene-annotation files (Human.GRCh38.p13.annot.tsv.gz etc.) are not served at
    # this URL (404 as of 2026-11-04); GeneIDs are mapped with org.Hs.eg.db / org.Mm.eg.db
    # in Scripts/process_cnp_source_data.R instead.
    got_counts <- FALSE
    for (f in gen) {
      is_annot <- grepl("annot", f)
      if (is_annot && !got_counts) next     # only take annotation when counts exist
      genome <- if (grepl("GRCh38", f)) "GRCh38" else "GRCm39"
      if (is_annot && !any(grepl(genome, list.files(file.path(dest_dir, "ncbi_counts"))))) next
      url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s", acc, f)
      ok <- fetch(url, file.path(dest_dir, "ncbi_counts", f), group, acc, role, quiet_404 = TRUE)
      if (isTRUE(ok) && !is_annot) got_counts <- TRUE
    }
    if (!got_counts) logmsg("  (no NCBI-generated counts for ", acc, ")")
  }
}

download_gds <- function(acc, group, role) {
  logmsg(acc, " (", group, ")")
  base <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/datasets/%s/%s/soft", geo_stub(acc), acc)
  for (f in c(paste0(acc, "_full.soft.gz"), paste0(acc, ".soft.gz")))
    fetch(file.path(base, f), file.path(OUT_DIR, group, acc, f), group, acc, role)
}

# ---- GWAS Catalog (REST API v2) ---------------------------------------------
GWAS_API <- "https://www.ebi.ac.uk/gwas/rest/api/v2"

collapse_lists <- function(df) {
  for (nm in names(df)) if (is.list(df[[nm]]))
    df[[nm]] <- vapply(df[[nm]], function(x) paste(unlist(x), collapse = ";"), "")
  df
}

gwas_get_all <- function(endpoint, params, max_pages = 500) {
  pages <- list(); page <- 0
  repeat {
    q <- paste0(names(params), "=", vapply(params, URLencode, "", reserved = TRUE), collapse = "&")
    url <- sprintf("%s/%s?%s&size=200&page=%d", GWAS_API, endpoint, q, page)
    js <- tryCatch(fromJSON(paste(readLines(url, warn = FALSE), collapse = "\n"), simplifyVector = FALSE),
                   error = function(e) { logmsg("  GWAS request failed: ", url, " (", conditionMessage(e), ")"); NULL })
    if (is.null(js)) break
    pages[[length(pages) + 1]] <- js
    total_pages <- if (!is.null(js$page$totalPages)) js$page$totalPages else 1
    page <- page + 1
    if (page >= total_pages || page >= max_pages) break
    Sys.sleep(0.15)                     # stay under the 15 requests/s limit
  }
  pages
}

save_gwas <- function(pages, embedded_key, dest_stub, group, accession, role, url_desc) {
  dir.create(dirname(dest_stub), recursive = TRUE, showWarnings = FALSE)
  json_file <- paste0(dest_stub, ".json")
  write_json(pages, json_file, auto_unbox = TRUE, pretty = TRUE, null = "null")
  # Pool records from every page first, so pages with differing fields still combine
  records <- unlist(lapply(pages, function(p) p[["_embedded"]][[embedded_key]]), recursive = FALSE)
  rows <- if (length(records))
    collapse_lists(fromJSON(toJSON(records, auto_unbox = TRUE, null = "null"), flatten = TRUE))
  else NULL
  n <- if (is.null(rows)) 0 else nrow(rows)
  if (n > 0) {
    if ("association_id" %in% names(rows)) rows <- rows[!duplicated(rows$association_id), ]
    write.table(rows, paste0(dest_stub, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    add_manifest(group, accession, paste0(dest_stub, ".tsv"), url_desc, sprintf("downloaded (%d rows)", nrow(rows)), role)
  } else {
    add_manifest(group, accession, json_file, url_desc, "no records returned - check JSON", role)
  }
  add_manifest(group, accession, json_file, url_desc, "raw API response", role)
  logmsg("  ", accession, ": ", n, " ", embedded_key)
}

# ---- Run --------------------------------------------------------------------
logmsg("==== CNP source-data download started; output: ", OUT_DIR)

for (i in seq_len(nrow(geo))) {
  acc <- geo$accession[i]
  if (acc %in% processed_acc) { logmsg(acc, ": already processed into data/processed - skipped"); next }
  if (acc %in% EXCLUDED) { logmsg(acc, ": excluded from validation - skipped (see EXCLUDED)"); next }
  tryCatch(
    if (startsWith(acc, "GDS")) download_gds(acc, geo$group[i], geo$role[i])
    else download_gse(acc, geo$group[i], geo$role[i]),
    error = function(e) logmsg("  ERROR on ", acc, ": ", conditionMessage(e)))
  save_manifest()
}

logmsg("---- GWAS Catalog studies")
for (i in seq_len(nrow(gwas_studies))) {
  acc <- gwas_studies$accession[i]; grp <- gwas_studies$group[i]; role <- gwas_studies$role[i]
  stub <- file.path(OUT_DIR, grp, "GWAS", acc)
  tryCatch({
    save_gwas(gwas_get_all("studies", list(accession_id = acc)), "studies",
              paste0(stub, "_study_metadata"), grp, acc, role, sprintf("%s/studies?accession_id=%s", GWAS_API, acc))
    save_gwas(gwas_get_all("associations", list(accession_id = acc)), "associations",
              paste0(stub, "_associations"), grp, acc, role, sprintf("%s/associations?accession_id=%s", GWAS_API, acc))
  }, error = function(e) {
    logmsg("  ERROR on ", acc, ": ", conditionMessage(e))
    add_manifest(grp, acc, stub, GWAS_API, paste("failed:", conditionMessage(e)), role)
  })
  save_manifest()
}

logmsg("---- GWAS Catalog traits")
for (i in seq_len(nrow(gwas_traits))) {
  acc <- gwas_traits$accession[i]; grp <- gwas_traits$group[i]; role <- gwas_traits$role[i]
  stub <- file.path(OUT_DIR, grp, "GWAS", paste0(acc, "_trait"))
  tryCatch(
    save_gwas(gwas_get_all("associations", list(efo_id = acc)), "associations",
              paste0(stub, "_associations"), grp, acc, role, sprintf("%s/associations?efo_id=%s", GWAS_API, acc)),
    error = function(e) {
      logmsg("  ERROR on ", acc, ": ", conditionMessage(e))
      add_manifest(grp, acc, stub, GWAS_API, paste("failed:", conditionMessage(e)), role)
    })
  save_manifest()
}

save_manifest()
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# ---- Summary ----------------------------------------------------------------
cat("\n==== Summary ====\n")
print(table(manifest$group, sub(" .*|:.*", "", manifest$status)))
problems <- manifest[grepl("failed|skipped|no records|not available", manifest$status), c("accession", "file", "status")]
if (nrow(problems)) { cat("\nNeeds attention:\n"); print(problems, row.names = FALSE) }
logmsg("==== Done. Manifest: ", manifest_file)
