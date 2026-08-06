############################################################
# Autoimmune × Diet Microbiome — End-to-End Pipeline (FULL)
# Top 4 BioProjects: IBDMDB, RA, MS-16S, Celiac (GFD)
# Author: Kiloni Quiles
# Last updated: 2026-03-16 (ET)
#
# WHAT IT DOES
#  1) Acquisition: SRA/ENA download (prefetch -> fasterq-dump) [shotgun/16S]
#  2) Profiling: MetaPhlAn 4 (shotgun) + optional HUMAnN 3 pathways
#  3) Profiling: QIIME 2 (DADA2) + SILVA taxonomy (16S)
#  4) Integration: build species/SGB tables, harmonize metadata (incl. diet tags)
#  5) Associations: MaAsLin2 (multivariable); ANCOM-BC2 (confirmatory DA)
#  6) Correlations: SparCC; Networks: SPIEC-EASI
#  7) ML: SIAMCAT cross-validated disease classifiers (+ ROC/PR curves)
#  8) (Optional) Batch correction: ComBat-Seq for count-like features
#  9) Manuscript figures + CSVs (300-dpi PNG + vector PDF)
#
# CITATIONS (methods/choices)
#  - SRA Run Selector + AccessionList/RunInfo best practice         (NCBI)         https://trace.ncbi.nlm.nih.gov/Traces/study1/?go=help
#  - SRA Toolkit prefetch/fasterq-dump usage and HPC notes           (NCBI/NIH HPC) https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump ; https://hpc.nih.gov/apps/sratoolkit.html
#  - MetaPhlAn 4 (species/SGB profiling)                            (bioBakery)     https://forum.biobakery.org/t/metaphlan-4-published-database-update/4850
#  - HUMAnN 3 (pathways; MetaCyc/UniRef)                            (bioBakery)     https://huttenhower.sph.harvard.edu/humann/
#  - QIIME 2 / DADA2 + SILVA taxonomy                               (QIIME2 docs)   https://docs.qiime2.org/2024.10/plugins/available/dada2/index.html
#  - MaAsLin 2 (multivariable associations)                          (Bioconductor)  https://bioconductor.org/packages//release/bioc/html/Maaslin2.html
#  - ANCOM-BC2 (bias-corrected DA)                                   (Bioconductor)  https://bioconductor.org/packages//release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
#  - SparCC (compositional correlation)                              (PLOS Comp Bio) https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687
#  - SPIEC-EASI (partial-correlations; networks)                     (SpiecEasi)     https://bioc.r-universe.dev/SpiecEasi/doc/manual.html
#  - SIAMCAT (ML toolbox for microbiomes)                            (Bioc)          https://bioconductor.posit.co/packages/3.19/bioc/html/SIAMCAT.html
#  - ComBat-Seq (NB batch correction on counts)                      (NAR G&B)       https://academic.oup.com/nargab/article/2/3/lqaa078/5909519
#  - Diet vocabularies: FoodOn; USDA WWEIA/FNDDS                     (FoodOn, USDA)  https://foodon.org/ ; https://www.ars.usda.gov/.../dmr-food-categories/
############################################################

## 0) USER CONFIGURATION ----
# Edit these paths and toggles to match your environment before running.

CONFIG <- list(
  
  # Root project folder (will be created if absent)
  ROOT = normalizePath("autoimmune_microbiome_fullrun", mustWork = FALSE),
  
  # Where to keep raw FASTQs (big!); prefer fast local scratch or RAID
  FASTQ_DIR = "~/scratch/autoimmune_fastq",
  
  # SRA Toolkit cache (vdb-config) — use a large scratch path if possible
  SRA_CACHE = "~/scratch/ncbi_sra_cache",
  
  # External tools (if installed on PATH, keep as-is)
  QIIME2_ENV = "~/miniconda3/envs/qiime2-amplicon",  # conda env name or prefix  # QIIME2/DADA2 (docs)  https://docs.qiime2.org/2024.10/...  [6](https://docs.qiime2.org/2024.10/plugins/available/dada2/index.html)
  METAPHLAN = "metaphlan",                            # MetaPhlAn 4 (bioBakery)   https://forum.biobakery.org/...                                    [4](https://www.mdpi.com/1424-8247/19/2/318)
  HUMANN    = "humann",                               # HUMAnN 3 (optional)        https://huttenhower.sph.harvard.edu/humann/                      [5](https://huttenhower.sph.harvard.edu/humann/)
  
  # Threading
  THREADS = 12,
  
  # Databases (only if running HUMAnN)
  HUMANN_DB = list(UNIREFFLAGS = "--protein-database uniref90", METACYC = "metacyc"),   # MetaCyc ref.  [27](https://metacyc.org/)
  
  # BioProjects (top 4)
  STUDIES = list(
    IBDMDB  = list(ID="PRJNA398089", TYPE="WGS",  STUDY="IBDMDB"),
    RA      = list(ID="PRJEB6997",   TYPE="WGS",  STUDY="RA"),     # ERP006678
    MS_16S  = list(ID="PRJNA321051", TYPE="16S",  STUDY="MS"),
    Celiac  = list(ID="PRJNA901394", TYPE="WGS",  STUDY="Celiac")
  ),
  
  # SILVA-trained classifier (QIIME2) paths (provide your own if you have them)
  SILVA_REF = list( # training details; see QIIME2/SILVA tutorial  [7](https://ejongepier.github.io/metabarcoding-qiime2-workshop/taxonomy.html)
    SEQS_QZA = "db/SILVA_138_99_16S-ref-seqs.qza",
    TAX_QZA  = "db/SILVA_138_99_16S-ref-taxonomy.qza"
  ),
  
  # Diet ontology mapping files (CSV you’ll curate; see FoodOn/WWEIA)
  DIET_MAPS = list(
    FOODON     = "metadata_maps/foodon_terms.csv",       # term, label, regex
    WWEIA_FNDDS= "metadata_maps/wweia_fndds_map.csv"     # code, category, regex
  ),
  
  # Optional: Use processed taxonomic tables for IBDMDB/MS to accelerate (set TRUE to skip raw re-profiling)
  USE_PROCESSED = list(IBDMDB=TRUE, MS_16S=TRUE)  # IBDMDB portal merged MetaPhlAn tables; MS paper supplement  [28](https://academic.oup.com/gastro/article/doi/10.1093/gastro/goae033/7659810)[29](https://www.researchgate.net/publication/345921954_ComBat-seq_batch_effect_adjustment_for_RNA-seq_count_data/fulltext/600c761392851c13fe32079a/ComBat-seq-batch-effect-adjustment-for-RNA-seq-count-data.pdf)
)

## 1) SETUP ----
dir.create(CONFIG$ROOT, recursive = TRUE, showWarnings = FALSE)
op <- options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  pkgs <- c(
    "data.table","dplyr","tidyr","readr","stringr","purrr",
    "ggplot2","pheatmap","vegan","igraph","forcats",
    # analysis
    "Maaslin2","ANCOMBC","SpiecEasi","SIAMCAT","phyloseq","sva"
  )
  BiocManager::install(pkgs[!pkgs %in% rownames(installed.packages())], ask = FALSE, update = TRUE)
  library(data.table); library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
  library(ggplot2); library(pheatmap); library(vegan); library(igraph); library(forcats)
  library(Maaslin2); library(ANCOMBC); library(SpiecEasi); library(SIAMCAT); library(phyloseq); library(sva)
})

# output folders
OUT <- file.path(CONFIG$ROOT, "outputs"); dir.create(OUT, showWarnings = FALSE)
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE)
CSV <- file.path(OUT, "tables");  dir.create(CSV, showWarnings = FALSE)
LOG <- file.path(CONFIG$ROOT, "logs"); dir.create(LOG, showWarnings = FALSE)
RPT <- file.path(OUT, "reports"); dir.create(RPT, showWarnings = FALSE)
dir.create(CONFIG$FASTQ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$SRA_CACHE, recursive = TRUE, showWarnings = FALSE)

# helper: save figures as PNG(300dpi) + PDF
save_figure <- function(p, base, w=10, h=6){
  ggsave(file.path(FIG, paste0(base,".png")), p, width=w, height=h, dpi=300, bg="white")
  ggsave(file.path(FIG, paste0(base,".pdf")), p, width=w, height=h, device=cairo_pdf)
}

## 2) DATA ACQUISITION (SRA/ENA) ----
# Best practice: use NCBI SRA Run Selector to export AccessionList and RunInfo TSV per BioProject,
# then feed that list to 'prefetch' (cache) and 'fasterq-dump' (FASTQ). [1](https://trace.ncbi.nlm.nih.gov/Traces/study1/?go=help)
# For HPC tips and tmp-space considerations see NIH HPC page and sra-tools wiki. [20](https://hpc.nih.gov/apps/sratoolkit.html)[2](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)

# Place your AccessionList files in 'accessions/' as {Study}_accessions.txt (one SRR/ERR per line)
ACC_DIR <- file.path(CONFIG$ROOT, "accessions"); dir.create(ACC_DIR, showWarnings = FALSE)

# Wrapper to download + convert an accession list to FASTQ
download_accessions <- function(acc_file, threads=CONFIG$THREADS){
  accs <- readLines(acc_file); accs <- accs[nchar(accs) > 0]
  if (length(accs) == 0) return(invisible(NULL))
  message("Prefetching ", length(accs), " runs... (", basename(acc_file), ")")
  # Use prefetch to cache .sra locally; more robust & resumable. [3](https://intro-to-bioinformatics-software.readthedocs.io/en/latest/data-acquisition.html)
  system2("prefetch", c("--output-directory", CONFIG$SRA_CACHE, "--option-file", acc_file), stdout=TRUE, stderr=TRUE)
  # Convert to FASTQ (paired split) with tmp on fast disk; see sra-tools wiki. [2](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
  for (acc in accs){
    cmd <- c(acc, "--split-files", "--threads", threads, "--progress", "--outdir", CONFIG$FASTQ_DIR,
             "--temp", CONFIG$SRA_CACHE)
    system2("fasterq-dump", cmd, stdout=TRUE, stderr=TRUE)
  }
}

# Example: run downloads if you’ve exported lists
# download_accessions(file.path(ACC_DIR,"IBDMDB_accessions.txt"))    # PRJNA398089   [22](https://www.cghjournal.org/article/S1542-3565%2824%2901048-6/fulltext)
# download_accessions(file.path(ACC_DIR,"RA_accessions.txt"))        # PRJEB6997/ERP006678 [23](https://www.cell.com/cell/fulltext/S0092-8674%2821%2900754-6)
# download_accessions(file.path(ACC_DIR,"MS_accessions.txt"))        # PRJNA321051   [25](https://www.biorxiv.org/content/10.1101/2020.01.13.904730v1.full.pdf)
# download_accessions(file.path(ACC_DIR,"Celiac_accessions.txt"))    # PRJNA901394   [26](https://link.springer.com/article/10.1186/s41927-025-00541-8)

## 3) PROFILING — Shotgun (MetaPhlAn 4) ----
# MetaPhlAn 4 yields species/SGB relative abundances (TSV). [4](https://www.mdpi.com/1424-8247/19/2/318)
# This block scans FASTQ_DIR and runs metaphlan per sample if profile missing.
run_metaphlan_for_dir <- function(fq_dir=CONFIG$FASTQ_DIR, pattern="_1.fastq.gz|_R1.fastq.gz"){
  fq1 <- list.files(fq_dir, pattern=pattern, full.names=TRUE)
  if (!length(fq1)) return(invisible(NULL))
  message("MetaPhlAn: profiling ", length(fq1), " samples...")
  out_dir <- file.path(OUT,"metaphlan"); dir.create(out_dir, showWarnings = FALSE)
  for (r1 in fq1){
    sample <- sub("_R1.*|_1.*","", basename(r1))
    r2 <- sub("_1.fastq.gz","_2.fastq.gz", r1)
    if (!file.exists(r2)) r2 <- sub("_R1.fastq.gz","_R2.fastq.gz", r1)
    profile <- file.path(out_dir, paste0(sample,"_profile.tsv"))
    if (file.exists(profile)) next
    args <- c(paste(r1, r2, sep=","), "--input_type","fastq", "--nproc", CONFIG$THREADS, "-o", profile)
    system2(CONFIG$METAPHLAN, args, stdout=TRUE, stderr=TRUE)
  }
  invisible(out_dir)
}

# (Optional) HUMAnN 3 functional profiling from the same FASTQs. [5](https://huttenhower.sph.harvard.edu/humann/)
run_humann_for_dir <- function(fq_dir=CONFIG$FASTQ_DIR){
  fq1 <- list.files(fq_dir, pattern="_1.fastq.gz|_R1.fastq.gz", full.names=TRUE)
  if (!length(fq1)) return(invisible(NULL))
  out_dir <- file.path(OUT,"humann"); dir.create(out_dir, showWarnings = FALSE)
  for (r1 in fq1){
    sample <- sub("_R1.*|_1.*","", basename(r1))
    r2 <- sub("_1.fastq.gz","_2.fastq.gz", r1); if (!file.exists(r2)) r2 <- sub("_R1.fastq.gz","_R2.fastq.gz", r1)
    if (!file.exists(file.path(out_dir, paste0(sample,"_genefamilies.tsv")))){
      args <- c("--input", r1, "--input2", r2, "--threads", CONFIG$THREADS, "--output", out_dir, CONFIG$HUMANN_DB$UNIREFFLAGS)
      system2(CONFIG$HUMANN, args, stdout=TRUE, stderr=TRUE)
    }
  }
  invisible(out_dir)
}

## 4) PROFILING — 16S (QIIME 2 / DADA2) ----
# This block *orchestrates* QIIME 2 steps via system calls from R (import, denoise, taxonomy using SILVA).
# See QIIME 2 DADA2 plugin docs. [6](https://docs.qiime2.org/2024.10/plugins/available/dada2/index.html)
qiime_run_16S <- function(manifest_csv, paired=TRUE){
  qiime <- function(...) system2(file.path(CONFIG$QIIME2_ENV,"bin","qiime"), c(...), stdout=TRUE, stderr=TRUE)
  qdir <- file.path(OUT,"qiime2_ms"); dir.create(qdir, showWarnings = FALSE)
  # 1) Import
  demux <- file.path(qdir,"demux.qza")
  if (!file.exists(demux)){
    type <- if (paired) "SampleData[PairedEndSequencesWithQuality]" else "SampleData[SequencesWithQuality]"
    fmt  <- if (paired) "PairedEndFastqManifestPhred33V2" else "SingleEndFastqManifestPhred33V2"
    qiime("tools","import","--type",type,"--input-path",manifest_csv,"--output-path",demux,"--input-format",fmt)
  }
  # 2) DADA2 denoise
  tbl <- file.path(qdir,"table.qza"); rep <- file.path(qdir,"rep_seqs.qza")
  if (!file.exists(tbl)){
    if (paired){
      qiime("dada2","denoise-paired","--i-demultiplexed-seqs",demux,"--p-trunc-len-f","0","--p-trunc-len-r","0",
            "--o-table",tbl,"--o-representative-sequences",rep,"--o-denoising-stats",file.path(qdir,"denoise_stats.qza"))
    } else {
      qiime("dada2","denoise-single","--i-demultiplexed-seqs",demux,"--p-trunc-len","0",
            "--o-table",tbl,"--o-representative-sequences",rep,"--o-denoising-stats",file.path(qdir,"denoise_stats.qza"))
    }
  }
  # 3) Taxonomy with SILVA
  taxa <- file.path(qdir,"taxonomy.qza")
  if (!file.exists(taxa)){
    qiime("feature-classifier","classify-sklearn","--i-classifier",CONFIG$SILVA_REF$TAX_QZA,"--i-reads",rep,"--o-classification",taxa)
  }
  # 4) Export to CSV
  system2(file.path(CONFIG$QIIME2_ENV,"bin","qiime"), c("tools","export","--input-path",tbl,"--output-path",qdir))
  system2(file.path(CONFIG$QIIME2_ENV,"bin","qiime"), c("tools","export","--input-path",taxa,"--output-path",qdir))
  invisible(qdir)
}

## 5) INTEGRATION — build abundance tables (species-level) + metadata ----
# Utilities to read MetaPhlAn tsvs -> merged table; read QIIME exports -> merged table.

read_metaphlan_merged <- function(path){
  files <- list.files(path, pattern="_profile.tsv$", full.names = TRUE); if (!length(files)) return(NULL)
  lst <- lapply(files, function(f){
    df <- fread(f, sep="\t"); # clade_name, relative_abundance, etc.
    df %>% filter(grepl("^s__", clade_name)) %>% select(Taxon=clade_name, Abundance=relative_abundance) %>%
      mutate(Sample = sub("_profile.tsv","", basename(f)))
  })
  long <- bind_rows(lst)
  wide <- long %>% pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
  as.data.frame(wide)
}

read_qiime_16S_species <- function(qiime_export_dir){
  # QIIME exported BIOM -> feature-table.biom; we assume table.tsv (feature x sample) also present
  tsv <- file.path(qiime_export_dir,"feature-table.tsv")
  tax <- file.path(qiime_export_dir,"taxonomy.tsv")
  if (!file.exists(tsv) || !file.exists(tax)) return(NULL)
  tbl <- fread(tsv, skip=1) %>% rename(FeatureID = `#OTU ID`)
  tx  <- fread(tax) %>% rename(FeatureID = `Feature ID`, Taxon = `Taxon`)
  dat <- tbl %>% inner_join(tx, by="FeatureID")
  # Collapse to species level (where available)
  dat$Species <- ifelse(grepl("s__", dat$Taxon), sub(".*;s__","s__", dat$Taxon), NA)
  dat <- dat %>% filter(!is.na(Species)) %>% select(-FeatureID, -Taxon) %>% rename(Taxon = Species)
  # Sum counts by species; then TSS-normalize to relative abundance
  abun <- dat %>% group_by(Taxon) %>% summarise(across(-Taxon, ~sum(.x, na.rm=TRUE)))
  samp_cols <- setdiff(colnames(abun),"Taxon")
  abun[samp_cols] <- sweep(abun[samp_cols], 2, colSums(abun[samp_cols, , drop=FALSE]), "/") * 100
  as.data.frame(abun)
}

# Harmonize metadata and diet tags using FoodOn/WWEIA mappings provided by the user. [18](https://foodon.org/)[19](https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/dmr-food-categories/)
harmonize_metadata <- function(meta_csv, foodon_map=CONFIG$DIET_MAPS$FOODON, wweia_map=CONFIG$DIET_MAPS$WWEIA_FNDDS){
  md <- fread(meta_csv)
  if (file.exists(foodon_map)){
    m <- fread(foodon_map) # columns: regex,label,term
    for (i in seq_len(nrow(m))){
      hit <- grepl(m$regex[i], md$diet_notes, ignore.case = TRUE)
      md[[paste0("FOODON_", m$term[i])]] <- as.integer(hit)
    }
  }
  if (file.exists(wweia_map)){
    w <- fread(wweia_map) # regex -> WWEIA category flags
    for (i in seq_len(nrow(w))){
      hit <- grepl(w$regex[i], md$diet_notes, ignore.case = TRUE)
      md[[paste0("WWEIA_", w$category[i])]] <- as.integer(hit)
    }
  }
  md
}

## 6) ASSOCIATIONS — MaAsLin2 & ANCOM‑BC2 ----
# Multivariable models (Study random intercept), FDR BH. [10](https://bioconductor.org/packages//release/bioc/html/Maaslin2.html)
run_maaslin2 <- function(abundance_csv, metadata_csv, outdir){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  Maaslin2(input_data = abundance_csv,
           input_metadata = metadata_csv,
           output = outdir,
           fixed_effects = c("Disease","Age","Sex","Antibiotics","Diet_GFD","Diet_Fiber","Diet_Fat"),  # edit to your metadata
           random_effects = c("Study"),
           normalization = "TSS", transform = "LOG", analysis_method = "LM",
           max_significance = 0.25, correction = "BH", cores = CONFIG$THREADS)
  invisible(outdir)
}

# Confirmatory DA using ANCOM‑BC2 (supports random effects and multi‑group). [11](https://bioconductor.org/packages//release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html)
run_ancombc2 <- function(ps, fix_formula = "Disease + Age + Sex + Antibiotics",
                         rand_formula="(1|Study)", out_csv){
  res <- ANCOMBC::ancombc2(data = ps,
                           fix_formula = fix_formula, rand_formula = rand_formula,
                           p_adj_method = "holm", prv_cut = 0.1, lib_cut = 0,
                           alpha = 0.1, n_cl = max(1, floor(CONFIG$THREADS/2)))
  # Extract result table
  da <- res$res %>% as.data.frame()
  fwrite(da, out_csv)
  invisible(da)
}

## 7) CORRELATIONS (SparCC) & NETWORKS (SPIEC‑EASI) ----
# SpiecEasi includes a SparCC wrapper; we filter to top‑prevalent taxa and compute. [13](https://bioc.r-universe.dev/SpiecEasi/doc/manual.html)
run_sparcc <- function(abun, min_prev=0.2, out_csv=file.path(CSV,"correlations_SparCC.csv")){
  # abun: Taxon x samples (relative abundances)
  mat <- abun %>% tibble::column_to_rownames("Taxon") %>% as.matrix()
  prev <- rowMeans(mat > 0) ; keep <- names(prev[prev >= min_prev])
  M <- mat[keep, , drop=FALSE]
  cor <- SpiecEasi::sparcc(t(M))$Cor # samples in rows; transpose back
  C <- reshape2::melt(cor, varnames=c("Taxon1","Taxon2"), value.name = "r")
  C <- C %>% filter(Taxon1 < Taxon2)
  fwrite(C, out_csv)
  invisible(C)
}

# SPIEC‑EASI partial‑correlation network on top‑prevalent taxa. [13](https://bioc.r-universe.dev/SpiecEasi/doc/manual.html)
run_spieceasi <- function(abun, min_prev=0.2, out_edges=file.path(CSV,"network_SPIECEASI_edges.csv")){
  mat <- abun %>% tibble::column_to_rownames("Taxon") %>% as.matrix()
  prev <- rowMeans(mat > 0) ; keep <- names(prev[prev >= min_prev])
  M <- mat[keep, , drop=FALSE]
  se <- spiec.easi(t(M), method='mb', lambda.min.ratio=1e-2, nlambda=20)
  B  <- as.matrix(symBeta(getOptBeta(se), mode='maxabs'))
  E  <- reshape2::melt(B, varnames=c("Taxon1","Taxon2"), value.name = "weight") %>%
    filter(Taxon1 < Taxon2, weight != 0)
  fwrite(E, out_edges)
  invisible(E)
}

## 8) MACHINE LEARNING — SIAMCAT (Disease classification) ----
# Build cross‑validated LASSO models across studies. [14](https://bioconductor.posit.co/packages/3.19/bioc/html/SIAMCAT.html)
run_siamcat <- function(abun, meta, label_col="Disease", out_prefix="siamcat"){
  feat <- abun %>% tibble::column_to_rownames("Taxon") %>% as.matrix()
  md   <- meta %>% tibble::column_to_rownames("SampleID")
  common <- intersect(colnames(feat), rownames(md))
  feat <- feat[, common, drop=FALSE]; md <- md[common, , drop=FALSE]
  sc <- siamcat(feat = feat, meta = md, label = label(md, label_col, positive_class = "Case"))
  sc <- filter.features(sc, cutoff=0.0) %>% normalize.features("log.unit")
  sc <- create.data.split(sc, num.folds=10, num.resample=5, stratify = TRUE, seed = 2026)
  sc <- train.model(sc, method="lasso") %>% make.predictions() %>% evaluate.predictions()
  png(file.path(FIG, paste0(out_prefix, "_model_eval.png")), width=1800, height=1200, res=300)
  print(model.evaluation.plot(sc))
  dev.off()
  saveRDS(sc, file.path(OUT, paste0(out_prefix,".rds")))
  invisible(sc)
}

## 9) (OPTIONAL) BATCH ADJUSTMENT — ComBat‑Seq on counts (HUMAnN gene families) ----
# Only for count‑like features (e.g., HUMAnN gene family counts prior to normalization).
# See ComBat‑Seq paper for caveats. [16](https://academic.oup.com/nargab/article/2/3/lqaa078/5909519)
combat_seq_adjust <- function(count_mat, batch_vec, group=NULL){
  adj <- ComBat_seq(counts = as.matrix(count_mat), batch = batch_vec, group = group, full_mod = TRUE)
  return(adj)
}

## 10) FIGURES (Profiling/PCoA/Associations/Networks) ----
plot_stacked_bars <- function(abun, meta, study, group_var="Disease", topN=20){
  dat <- abun %>% pivot_longer(-Taxon, names_to="SampleID", values_to="Abundance") %>%
    inner_join(meta, by="SampleID") %>% filter(Study==study)
  keep <- dat %>% group_by(Taxon) %>% summarize(m=mean(Abundance), .groups="drop") %>% arrange(desc(m)) %>%
    slice_head(n=topN) %>% pull(Taxon)
  dat2 <- dat %>% filter(Taxon %in% keep) %>% group_by(SampleID) %>%
    mutate(rel=Abundance/sum(Abundance)) %>% ungroup()
  p <- ggplot(dat2, aes(x=SampleID, y=rel, fill=Taxon)) + geom_col() +
    facet_grid(~ .data[[group_var]], scales="free_x", space="free_x") +
    scale_y_continuous(labels=scales::percent_format()) +
    labs(x="Samples", y="Relative abundance", title=paste("Top taxa —", study)) +
    theme_classic(base_size=11) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  save_figure(p, paste0("stackedbar_",study))
}

plot_pcoa <- function(abun, meta){
  mat <- abun %>% tibble::column_to_rownames("Taxon") %>% as.matrix()
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  bc  <- vegan::vegdist(mat, method="bray")  # PCoA (Bray–Curtis)  [30](https://rstudio-pubs-static.s3.amazonaws.com/794628_9729e7e28b0049ab98b61adedbdda2a5.html)
  pc  <- cmdscale(bc, k=2)
  ord <- data.frame(pc) %>% tibble::rownames_to_column("SampleID") %>% inner_join(meta, by="SampleID")
  p <- ggplot(ord, aes(x=X1, y=X2, color=Disease, shape=Study)) + geom_point(size=2, alpha=0.9) +
    labs(x="PCoA1", y="PCoA2", title="PCoA (Bray–Curtis)") + theme_classic(base_size=11)
  save_figure(p, "pcoa_bray")
}

plot_assoc_heatmap <- function(maslin_results_tsv, title="MaAsLin2 β (q≤0.1)"){
  res <- fread(maslin_results_tsv)  # standard MaAsLin2 output
  top <- res %>% filter(qval <= 0.1) %>% select(feature, metadata, coef) %>%
    pivot_wider(names_from=metadata, values_from=coef, values_fill=0) %>%
    tibble::column_to_rownames("feature") %>% as.matrix()
  pheatmap(top, color=colorRampPalette(c("#3B4CC0","white","#B40426"))(101), main=title)
  ggsave(file.path(FIG,"assoc_maaslin_heatmap.pdf"), width=8.5, height=6, device=cairo_pdf)
  ggsave(file.path(FIG,"assoc_maaslin_heatmap.png"), width=8.5, height=6, dpi=300, bg="white")
}

plot_network <- function(edges_csv){
  E <- fread(edges_csv)
  g <- graph_from_data_frame(E, directed=FALSE)
  deg <- degree(g)
  V(g)$size <- scales::rescale(deg, to=c(3,12))
  V(g)$label <- NA
  E(g)$color <- ifelse(E(g)$weight >= 0, "#2C7BB6", "#D7191C")
  png(file.path(FIG,"network_spieceasi.png"), width=2200, height=1800, res=300)
  plot(g, layout=layout_with_fr(g), main="SPIEC‑EASI (node size = degree)")
  dev.off()
  pdf(file.path(FIG,"network_spieceasi.pdf"), width=10, height=8)
  plot(g, layout=layout_with_fr(g), main="SPIEC‑EASI (node size = degree)")
  dev.off()
}

## 11) PIPELINE DRIVER ----
# This section shows how to wire the pieces together end‑to‑end.
# Adjust 'meta_master.csv' to your metadata schema; required columns:
#  SampleID, Study, Disease, Age, Sex, Antibiotics, diet_notes, (optional diet harmonized flags)

# (A) Build species tables
build_abundance_tables <- function(){
  # Shotgun — either use processed IBDMDB or run MetaPhlAn
  mp_dir <- if (CONFIG$USE_PROCESSED$IBDMDB) {
    # drop your downloaded IBDMDB merged MetaPhlAn table here (from portal)  [28](https://academic.oup.com/gastro/article/doi/10.1093/gastro/goae033/7659810)
    # Otherwise call: run_metaphlan_for_dir()
    message("Using IBDMDB processed MetaPhlAn merged tables per portal.")
    # Example: read a prepared metaphlan folder:
    file.path(OUT,"metaphlan")  # assume ready
  } else run_metaphlan_for_dir()
  
  abun_shotgun <- if (dir.exists(file.path(OUT,"metaphlan"))) read_metaphlan_merged(file.path(OUT,"metaphlan")) else NULL
  
  # MS 16S — either use processed supplement or run QIIME 2
  if (CONFIG$USE_PROCESSED$MS_16S){
    message("Using MS processed 16S supplement for species table.")  # Jangi et al. 2016  [29](https://www.researchgate.net/publication/345921954_ComBat-seq_batch_effect_adjustment_for_RNA-seq_count_data/fulltext/600c761392851c13fe32079a/ComBat-seq-batch-effect-adjustment-for-RNA-seq-count-data.pdf)
    # Place exported table as CSV: Taxon + samples
    ms_abun <- fread(file.path(CONFIG$ROOT,"processed/ms16s_species.csv"))
  } else {
    qiime_run_16S(manifest_csv = file.path(CONFIG$ROOT,"ms_manifest.csv"), paired=TRUE)  # you provide a manifest
    ms_abun <- read_qiime_16S_species(file.path(OUT,"qiime2_ms"))
  }
  
  # Merge shotgun + ms where possible (outer join on Taxon)
  all_abun <- list()
  if (!is.null(abun_shotgun)) all_abun[["shotgun"]] <- abun_shotgun
  if (exists("ms_abun")) all_abun[["ms"]] <- as.data.frame(ms_abun)
  abun_merged <- reduce(all_abun, ~full_join(.x, .y, by="Taxon")) %>% replace(is.na(.), 0)
  
  # Persist
  fwrite(abun_merged, file.path(CSV,"abundance_taxa.csv"))
  return(abun_merged)
}

# (B) Run stats / figs
run_full_analytics <- function(){
  abun <- fread(file.path(CSV,"abundance_taxa.csv"))
  meta <- harmonize_metadata(file.path(CONFIG$ROOT,"meta_master.csv"))
  fwrite(meta, file.path(CSV,"sample-level_metadata.csv"))
  
  # Per-study stacked bars + merged PCoA
  for (s in c("IBDMDB","RA","Celiac","MS")){
    try(plot_stacked_bars(abun, meta, s), silent=TRUE)
  }
  try(plot_pcoa(abun, meta), silent=TRUE)
  
  # MaAsLin2 (merged)
  mas_dir <- file.path(OUT,"maaslin2"); run_maaslin2(file.path(CSV,"abundance_taxa.csv"), file.path(CSV,"sample-level_metadata.csv"), mas_dir)
  # Heatmap from MaAsLin2 all results (edit path as needed)
  res_tsv <- file.path(mas_dir,"all_results.tsv"); if (file.exists(res_tsv)) plot_assoc_heatmap(res_tsv)
  
  # ANCOM‑BC2 (build phyloseq from relative abundances for convenience)
  mat <- abun %>% tibble::column_to_rownames("Taxon") %>% as.matrix()
  ps  <- phyloseq(otu_table(mat, taxa_are_rows = TRUE), sample_data(meta %>% tibble::column_to_rownames("SampleID")))
  run_ancombc2(ps, out_csv = file.path(CSV,"ancombc2_disease_DA.csv"))
  
  # SparCC + SPIEC‑EASI
  run_sparcc(abun, min_prev=0.2, out_csv=file.path(CSV,"correlations_SparCC.csv"))
  run_spieceasi(abun, min_prev=0.2, out_edges=file.path(CSV,"network_SPIECEASI_edges.csv"))
  plot_network(file.path(CSV,"network_SPIECEASI_edges.csv"))
  
  # SIAMCAT ML
  run_siamcat(abun, meta, label_col="Disease", out_prefix="siamcat_species")
  
  message("All analytics complete. Tables in ", CSV, ", figures in ", FIG, ".")
}

# ------------- EXECUTE -------------
# 1) Build abundance tables (or re-run to refresh)
# abun_merged <- build_abundance_tables()

# 2) Run analytics and figure generation
# run_full_analytics()

############################################################
# END OF PIPELINE
############################################################