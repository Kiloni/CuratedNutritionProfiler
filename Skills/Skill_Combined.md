---
name: nutrient-signature-profiler

description: Discovery, harmonization, independent validation, host regulatory mapping, and precision nutrition interpretation of host–microbiome–diet multi-omics datasets. Use when analyzing autoimmune microbiome studies, validating microbial signatures, integrating host transcriptomics and microbiome data, performing cross-kingdom analyses, generating precision nutrition recommendations, or reproducing NutrientSignatureProfiler workflows.
---

# NutrientSignatureProfiler

A reusable workflow for discovery, harmonization, independent validation, host regulatory mapping, and precision nutrition interpretation of host–microbiome–diet datasets.

The workflow focuses on autoimmune disease microbiome research using public multi-omics resources and emphasizes validation before biomarker reporting.

---

## When to Trigger

Use this skill when:

- Building microbiome biomarker catalogues
- Validating autoimmune microbiome findings
- Integrating microbiome and transcriptomic datasets
- Analyzing BugSigDB signatures
- Analyzing GMrepo cohort fingerprints
- Performing host–microbe interaction analyses
- Reproducing NutrientSignatureProfiler workflows
- Generating precision nutrition recommendations
- Constructing validation-focused manuscripts

Primary disease targets:

- Multiple Sclerosis (MS)
- Systemic Lupus Erythematosus (SLE)
- Rheumatoid Arthritis (RA)
- Chronic Fatigue Syndrome (CFS)

---

## Core Principles

1. Preserve Original Project ID.
2. Separate discovery and validation cohorts.
3. Never report a signature as validated without independent replication.
4. Remove metabolic and technical confounders.
5. Preserve study provenance.
6. Save checkpoints after each phase.
7. Record filtering decisions.
8. Prioritize reproducibility over novelty.
9. Retain all validation evidence.
10. Generate publication-ready outputs.

---

## Scope

This skill provides workflows for:

- Cohort discovery
- Microbiome harmonization
- Host–microbiome integration
- Validation studies
- Pathway-level interpretation
- Precision nutrition recommendation generation

This skill is intended for research workflows and does not replace clinical decision making.

---

## Runtime Requirements

### Recommended R Environment

Reference implementation:

```r
install.packages(c(
  "tidyverse",
  "data.table",
  "janitor",
  "arrow",
  "plotly",
  "cowplot",
  "patchwork",
  "vegan",
  "randomForest",
  "shiny",
  "DT"
))

BiocManager::install(c(
  "BugSigDB",
  "MultiAssayExperiment",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "limma",
  "GSVA",
  "fgsea",
  "clusterProfiler",
  "ReactomePA",
  "GEOquery",
  "viper"
))
```

### Recommended External Tools

Reference implementation:

```bash
conda install -c bioconda sra-tools -y
conda install -c bioconda bwa samtools -y
conda install -c bioconda salmon -y
conda install -c bioconda entrez-direct -y
```

---

## Repository Structure

```text
NutrientSignatureProfiler/

├── data/
│   ├── host/
│   ├── microbiome/
│   ├── metadata/
│   └── validation/
│
├── outputs/
│   ├── phase1/
│   ├── phase2/
│   ├── phase3/
│   ├── phase4_validation/
│   ├── phase5_fidelity/
│   ├── phase6_nutrition/
│   └── reports/
│
├── scripts/
├── notebooks/
├── validation/
└── figures/
```

---

## Inputs

### Microbiome

- Abundance matrices
- LEfSe outputs
- GMrepo signatures
- BugSigDB signatures
- Cohort metadata
- Taxonomic profiles
- Functional profiles

### Host Omics

- Bulk RNA-seq
- Normalized expression matrices
- Differential expression results
- Celligner outputs
- Gene-level statistics

### Metadata

- Disease labels
- Cohort identifiers
- Subject identifiers
- Original Project IDs
- Technical metadata
- Geographic metadata

### Optional

- Dietary metadata
- Metabolomics
- Proteomics
- ARACNe networks
- metaVIPER regulons
- Clinical laboratory values

---

# Phase I — Discovery

## Purpose

Construct matched host–microbiome cohort inventories suitable for downstream harmonization and independent validation.

This phase is responsible for identifying:

- Discovery cohorts
- Validation cohorts
- Host transcriptomic datasets
- Microbiome cohorts
- Supporting metadata resources

All studies should maintain traceability through Original Project IDs.

Unlike downstream phases, this phase does not require uploaded matrices.

The focus is public cohort discovery, metadata retrieval, and validation dataset identification.

---

## Data Sources

Search:

- GEO
- SRA
- ENA
- GMrepo
- BugSigDB
- CuratedMetagenomicData

Do not hardcode accession identifiers.

Search repositories dynamically so future datasets can automatically be incorporated into the workflow.

---

## Workflow

1. Search GEO.
2. Search SRA.
3. Search ENA.
4. Search GMrepo.
5. Search BugSigDB.
6. Retrieve real study metadata.
7. Build candidate catalogues.
8. Identify validation-eligible cohorts.
9. Generate discovery metadata inventories.
10. Export validation catalogue.

Do not fabricate:

- GEO accessions
- SRA accessions
- Sample counts
- Study metadata
- Platform information

---

## Discovery Rules

Primary disease targets:

- Multiple Sclerosis
- Systemic Lupus Erythematosus
- Rheumatoid Arthritis
- Chronic Fatigue Syndrome

Secondary disease targets:

- Additional autoimmune disorders
- Immune-mediated chronic diseases

---

## Required Metadata Fields

Capture whenever available:

- Study Accession
- Original Project ID
- Sample ID
- Subject ID
- Disease
- Tissue Source
- Microbiome Source
- Platform
- Country
- Read Length
- Sample Count
- PMID

---

## Comorbidity Schema

| Column | Description |
|----------|-------------|
| A | Cancer Comorbidities |
| B | Metabolic Comorbidities |
| C | Autoimmune Comorbidities |
| D | Other Comorbidities |

---

## Validation Dataset Discovery Rules

Validation cohorts should:

- Be independent from discovery cohorts
- Use human subjects
- Contain microbiome profiles
- Preferably contain host transcriptomics
- Preferably contain dietary metadata
- Have adequate sample counts
- Preserve accession traceability

Validation cohorts must not be randomly generated.

Validation eligibility should be based on retrieved study metadata.

---

## Reference Implementation

```r
############################################################
# Phase I — Public Cohort Discovery
############################################################

library(GEOquery)
library(dplyr)
library(stringr)
library(tibble)

############################################################
# Disease Search Terms
############################################################

disease_queries <- c(

  "multiple sclerosis microbiome",

  "systemic lupus erythematosus microbiome",

  "rheumatoid arthritis microbiome",

  "chronic fatigue syndrome microbiome",

  "ME/CFS microbiome"

)

############################################################
# Candidate Study Catalogue
############################################################

geo_candidate_catalog <- list()

############################################################
# Discovery Notes
############################################################

message(
  "Beginning Phase I discovery."
)

message(
  "Search GEO, SRA, ENA, GMrepo and BugSigDB."
)

message(
  "Do not create placeholder accession IDs."
)

message(
  "Only record studies that actually exist."
)

############################################################
# GEO Retrieval Template
############################################################

#
# Example workflow:
#
# 1. Search GEO manually or through repository APIs.
# 2. Identify a valid GSE accession.
# 3. Retrieve metadata.
# 4. Append study information.
#
# Example:
#
# gse <- getGEO("GSEXXXXXX")
#
# metadata <- pData(
#   phenoData(
#     gse[[1]]
#   )
# )
#
# study_row <- data.frame(
#
#     accession =
#       "GSEXXXXXX",
#
#     title =
#       Meta(gse)$title,
#
#     platform =
#       Meta(gse)$platform_id,
#
#     sample_count =
#       nrow(metadata),
#
#     stringsAsFactors = FALSE
#
# )
#
# geo_candidate_catalog[[length(
#   geo_candidate_catalog
# ) + 1]] <- study_row
#

############################################################
# Validation Eligibility Function
############################################################

validation_filter <- function(
  sample_count,
  microbiome_available,
  organism = "Homo sapiens"
){

  pass <- TRUE

  if(
    organism !=
    "Homo sapiens"
  ){

    pass <- FALSE

  }

  if(
    microbiome_available == FALSE
  ){

    pass <- FALSE

  }

  if(
    sample_count < 20
  ){

    pass <- FALSE

  }

  return(pass)

}

############################################################
# Candidate Catalogue Assembly
############################################################

if(
  length(
    geo_candidate_catalog
  ) > 0
){

  geo_tbl <- bind_rows(
    geo_candidate_catalog
  )

}else{

  geo_tbl <- tibble(

    accession = character(),

    title = character(),

    platform = character(),

    sample_count = numeric()

  )

}

############################################################
# Add Metadata Flags
############################################################

if(
  nrow(geo_tbl) > 0
){

  geo_tbl <- geo_tbl %>%

    mutate(

      validation_eligible =

        mapply(

          validation_filter,

          sample_count =
            sample_count,

          microbiome_available =
            TRUE

        ),

      discovery_status =
        "Discovered"

    )

}

############################################################
# Export Catalogues
############################################################

write.csv(

  geo_tbl,

  "Phase_I_GEO_Candidate_Studies.csv",

  row.names = FALSE

)

validation_tbl <-

  geo_tbl %>%

  filter(
    validation_eligible
  )

write.csv(

  validation_tbl,

  "Phase_I_Validation_Catalogue.csv",

  row.names = FALSE

)

############################################################
# Discovery Summary
############################################################

summary_tbl <- data.frame(

  total_candidates =
    nrow(geo_tbl),

  validation_candidates =
    nrow(validation_tbl),

  generated =
    Sys.time()

)

write.csv(

  summary_tbl,

  "Phase_I_Discovery_Summary.csv",

  row.names = FALSE

)

message(
  "Phase I complete."
)

message(
  "No placeholder accessions were generated."
)
```

---

## Execution Requirements

- Search public repositories before creating outputs.
- Do not generate fictional GEO accessions.
- Do not generate fictional SRA projects.
- Do not create NA-filled study catalogs.
- Only include studies that were actually retrieved.
- If no qualifying studies are found, return an empty catalogue and report that no studies met the criteria.
- Discovery cohorts and validation cohorts must be generated from real metadata.
- Validation eligibility must be determined from retrieved study information.

---

## Standard Outputs

```text
Phase_I_GEO_Candidate_Studies.csv

Phase_I_SRA_Candidate_Studies.csv

Phase_I_Validation_Catalogue.csv

Phase_I_Discovery_Summary.csv

Phase_I_Host_Metadata.csv

Phase_I_Microbiome_Metadata.csv
```

---

# Phase II — Harmonization & Metabolic Noise Isolation

## Purpose

Generate harmonized host and microbiome signatures while minimizing:

- Batch effects
- Platform effects
- Technical artifacts
- Metabolic confounding
- Diet-associated noise

This phase establishes disease-level signatures for downstream validation.

---

## Workflow

1. Retrieve raw host count matrices and microbiome abundance matrices, with subject-level covariates, from GEO/SRA.
2. Normalize host expression.
3. Correct batch effects.
4. Harmonize microbiome signatures.
5. Build weighted meta-signatures.
6. Remove metabolic confounders.
7. Construct cleaned signature objects.
8. Flag dietary-stress cluster taxa.
9. Export harmonized matrices.

---

## Raw Matrix Retrieval from GEO Supplementary Files

Purpose:

Retrieve the raw host count matrices, raw microbiome abundance matrices, and
subject-level covariates (age, sex, BMI, diet) that every later step in this
phase depends on, keyed off the accessions discovered in Phase I. This step
did not previously exist in this skill; every downstream step silently
assumed its output (`host_counts.csv`, a `metadata` object) already existed.

Reference implementation:

```r
library(GEOquery)
library(Biobase)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

geo_candidates <- read_csv("Phase_I_GEO_Candidate_Studies.csv")
sra_candidates <- read_csv("Phase_I_SRA_Candidate_Studies.csv")

fetch_geo_eset <- function(accession) {
  tryCatch(
    getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)[[1]],
    error = function(e) { message(accession, " -> FAILED: ", conditionMessage(e)); NULL }
  )
}

esets <- geo_candidates %>%
  pull(accession) %>%
  unique() %>%
  set_names() %>%
  map(fetch_geo_eset) %>%
  compact()

# Platform matters downstream: not every host accession is RNA-seq
# (e.g. NanoString panels report already-processed counts, not raw
# library-size-scaled reads). Tag it here so normalization can branch
# correctly instead of blindly applying an RNA-seq-only method.
platform_lookup <- geo_candidates %>%
  select(accession, platform)

is_rnaseq_platform <- function(x) str_detect(str_to_lower(x), "rna-?seq|hiseq|novaseq|nextseq")

build_host_long <- function(accession) {
  eset <- esets[[accession]]
  if (is.null(eset) || nrow(exprs(eset)) == 0) return(NULL)
  platform <- platform_lookup %>% filter(accession == !!accession) %>% pull(platform)
  as.data.frame(exprs(eset)) %>%
    tibble::rownames_to_column("feature_id") %>%
    pivot_longer(-feature_id, names_to = "sample_id", values_to = "value") %>%
    mutate(accession = accession, platform = platform, is_rnaseq = is_rnaseq_platform(platform))
}

host_long <- geo_candidates %>%
  pull(accession) %>% unique() %>%
  map(build_host_long) %>% compact() %>% bind_rows() %>%
  mutate(sample_id = paste(accession, sample_id, sep = "__"))  # avoid cross-study ID collisions

host_counts_wide <- host_long %>%
  select(feature_id, sample_id, value) %>%
  pivot_wider(names_from = sample_id, values_from = value) %>%
  tibble::column_to_rownames("feature_id")

write.csv(host_counts_wide, "host_counts.csv")

host_sample_platform_map <- host_long %>%
  distinct(sample_id, accession, platform, is_rnaseq)

fetch_geo_covariates <- function(accession) {
  eset <- esets[[accession]]
  if (is.null(eset)) return(NULL)
  pdata <- pData(eset)
  covariate_cols <- grep("age|sex|gender|bmi|body.?fat|diet", colnames(pdata), ignore.case = TRUE, value = TRUE)
  pdata %>%
    select(geo_accession, any_of(covariate_cols)) %>%
    mutate(accession = accession, sample_id = paste(accession, geo_accession, sep = "__"), study_id = accession)
}

phase2_sample_metadata <- geo_candidates %>%
  pull(accession) %>% unique() %>%
  map_dfr(fetch_geo_covariates)

write_csv(phase2_sample_metadata, "Phase_II_Sample_Metadata.csv")
```

Microbiome raw abundance retrieval (`microbe_counts.csv`) is not a single
file download the way host supplementary files are: GEO/SRA/ENA do not
standardize a processed abundance-table format the way array platforms do
for expression data. Resolve it either from the paper's own Supplementary
Data (check first) or by running a full 16S (DADA2/QIIME2) or shotgun
(MetaPhlAn4/HUMAnN3) pipeline against FASTQ files resolved from the linked
BioProject:

```r
library(rentrez)

sra_candidates <- read_csv("Phase_I_SRA_Candidate_Studies.csv")

fetch_bioproject_run_list <- function(bioproject_accession) {
  tryCatch(
    entrez_search(db = "sra", term = paste0(bioproject_accession, "[BioProject]"), retmax = 10000),
    error = function(e) { message(bioproject_accession, " -> FAILED: ", conditionMessage(e)); NULL }
  )
}

sra_run_lists <- sra_candidates %>%
  pull(accession) %>% unique() %>% set_names() %>%
  map(fetch_bioproject_run_list) %>% compact()
```

Standard output:

```text
host_counts.csv
Phase_II_Sample_Metadata.csv
```

---

## Host TPM Normalization

Purpose:

Generate normalized host transcriptomic matrices suitable for integrative
analyses. Only true RNA-seq samples are routed through this method — a
NanoString or microarray sample forced through `edgeR::cpm()` produces a
number, but not a methodologically valid one.

Reference implementation:

```r
library(edgeR)

counts <- read.csv(
  "host_counts.csv",
  row.names = 1
)

rnaseq_samples <- host_sample_platform_map %>% filter(is_rnaseq) %>% pull(sample_id)

if (length(rnaseq_samples) > 0) {

  dge <- DGEList(
    counts[, rnaseq_samples, drop = FALSE]
  )

  tpm_matrix <- cpm(
    dge,
    normalized.lib.sizes = TRUE
  )

} else {

  message("No RNA-seq-platform samples available -- non-RNA-seq samples ",
          "in host_counts.csv need a platform-appropriate normalization ",
          "(not edgeR/CPM) before they can be added here.")
  tpm_matrix <- counts[, character(0), drop = FALSE]

}

write.csv(tpm_matrix, "Phase_II_Host_Harmonized.csv")
```

Standard output:

```text
Phase_II_Host_Harmonized.csv
```

---

## Batch Effect Correction

Purpose:

Remove study-level technical variation prior to validation.

Reference implementation:

```r
library(limma)

batch_vector <- phase2_sample_metadata %>%
  filter(sample_id %in% colnames(tpm_matrix)) %>%
  arrange(match(sample_id, colnames(tpm_matrix))) %>%
  pull(study_id)

stopifnot(length(batch_vector) == ncol(tpm_matrix))

expr_corrected <- removeBatchEffect(
  tpm_matrix,
  batch = batch_vector
)

write.csv(expr_corrected, "Phase_II_Host_Harmonized.csv")
```

Standard output:

```text
Phase_II_Host_Harmonized.csv
```

---

## Harmonize Microbiome Signatures

Purpose:

Normalize compositional microbiome abundance data and remove study-level
batch effects, mirroring the host TPM-normalization-plus-batch-correction
pair above. This step previously had no reference implementation despite
being named in the Workflow list and having its own Standard Output —
regression was run directly on raw counts, which skips normalization
entirely.

Reference implementation:

```r
microbe_counts_path <- "microbe_counts.csv"

if (file.exists(microbe_counts_path)) {

  microbe_raw <- read.csv(microbe_counts_path, row.names = 1)  # taxa x samples

  # Total-sum-scaling to relative abundance, then a centered log-ratio
  # (CLR) transform with a small pseudocount for zeros -- standard
  # handling for compositional microbiome data.
  relabund <- sweep(microbe_raw, 2, colSums(microbe_raw), "/")
  pseudo <- relabund + 1e-6
  clr_matrix <- apply(pseudo, 2, function(x) log(x) - mean(log(x)))
  rownames(clr_matrix) <- rownames(microbe_raw)

  microbe_batch_vector <- phase2_sample_metadata %>%
    filter(sample_id %in% colnames(clr_matrix)) %>%
    arrange(match(sample_id, colnames(clr_matrix))) %>%
    pull(study_id)

  if (length(unique(microbe_batch_vector)) > 1) {
    microbe_harmonized <- removeBatchEffect(clr_matrix, batch = microbe_batch_vector)
  } else {
    microbe_harmonized <- clr_matrix  # nothing to batch-correct against
  }

  write.csv(microbe_harmonized, "Phase_II_Microbe_Harmonized.csv")

} else {

  message(microbe_counts_path, " does not exist -- see Raw Matrix ",
          "Retrieval section above. Microbiome harmonization, metabolic ",
          "noise regression, and dietary-stress cluster detection all ",
          "depend on this file and are skipped without it.")

}
```

Standard output:

```text
Phase_II_Microbe_Harmonized.csv
```

---

## Metabolic Noise Regression

Purpose:

Remove microbiome signals associated with BMI, body fat, and dietary fat
before biomarker construction. Runs on the harmonized matrix from the
previous step, not on raw counts.

Reference implementation:

```r
microbe_matrix <- t(microbe_harmonized)  # samples x taxa

microbe_metadata <- phase2_sample_metadata %>%
  filter(sample_id %in% rownames(microbe_matrix)) %>%
  arrange(match(sample_id, rownames(microbe_matrix)))

required_cols <- c("bmi", "body_fat", "dietary_fat", "age", "sex")
missing_cols <- setdiff(required_cols, colnames(microbe_metadata))
if (length(missing_cols) > 0) {
  stop("microbe_metadata is missing required covariate column(s): ",
       paste(missing_cols, collapse = ", "))
}

microbe_adjusted <- list()

for(feature in colnames(microbe_matrix)){

  fit <- lm(
    microbe_matrix[,feature] ~
      microbe_metadata$bmi +
      microbe_metadata$body_fat +
      microbe_metadata$dietary_fat +
      microbe_metadata$age +
      microbe_metadata$sex
  )

  microbe_adjusted[[feature]] <- resid(fit)

}

write.csv(as.data.frame(microbe_adjusted), "Phase_II_Cleaned_Signatures.csv")
```

Standard output:

```text
Phase_II_Cleaned_Signatures.csv
```

---

## Dietary-Stress Cluster Detection

Purpose:

Monitor taxa known to track broad metabolic/dietary state rather than
disease-specific biology, so a signature isn't reported as disease-specific
without checking it against these confounders first.

Monitor:

- Akkermansia muciniphila
- Bacteroides fragilis
- Bilophila wadsworthia

Assess associations with:

- Host adiposity
- Dietary fat
- Chronic fatigue
- Metabolic states

These organisms help distinguish disease biology from broad metabolic effects.

Reference implementation:

```r
dietary_stress_taxa <- c(
  "Akkermansia muciniphila",
  "Bacteroides fragilis",
  "Bilophila wadsworthia"
)

present_taxa <- intersect(dietary_stress_taxa, rownames(microbe_harmonized))

test_association <- function(taxon_abundance, covariate) {
  if (is.null(covariate) || length(unique(na.omit(covariate))) < 2) return(NA_real_)
  tryCatch(cor.test(taxon_abundance, as.numeric(covariate))$p.value, error = function(e) NA_real_)
}

dsc_metadata <- phase2_sample_metadata %>%
  filter(sample_id %in% colnames(microbe_harmonized)) %>%
  arrange(match(sample_id, colnames(microbe_harmonized)))

dietary_stress_flags <- purrr::map_dfr(present_taxa, function(taxon) {
  taxon_abundance <- as.numeric(microbe_harmonized[taxon, ])
  tibble::tibble(
    taxon = taxon,
    present = TRUE,
    assoc_bmi_p = test_association(taxon_abundance, dsc_metadata$bmi),
    assoc_dietary_fat_p = test_association(taxon_abundance, dsc_metadata$dietary_fat),
    assoc_chronic_fatigue_p = test_association(taxon_abundance, dsc_metadata$chronic_fatigue)
  )
})

missing_taxa <- setdiff(dietary_stress_taxa, present_taxa)
if (length(missing_taxa) > 0) {
  dietary_stress_flags <- dplyr::bind_rows(
    dietary_stress_flags,
    tibble::tibble(taxon = missing_taxa, present = FALSE,
                    assoc_bmi_p = NA_real_, assoc_dietary_fat_p = NA_real_,
                    assoc_chronic_fatigue_p = NA_real_)
  )
}

write_csv(dietary_stress_flags, "Phase_II_Dietary_Stress_Cluster_Flags.csv")
```

Standard output:

```text
Phase_II_Dietary_Stress_Cluster_Flags.csv
```

---

## BugSigDB Consensus Signature Construction

Purpose:

Construct weighted disease meta-signatures from per-study log-odds-ratios.
Requires a `signatures` table with real `lor` and `sample_size` columns —
e.g. pulled from the BugSigDB API — which is a different object from a
literature-direction-only table (taxon/disease/direction, no `lor` or
`sample_size`) that a discovery-stage literature search alone produces.

Reference implementation:

```r
library(dplyr)
library(readr)

bugsigdb_signatures_path <- "Phase_II_BugSigDB_Signatures.csv"

if (file.exists(bugsigdb_signatures_path)) {

  signatures <- read_csv(bugsigdb_signatures_path)

  required_cols <- c("taxon", "disease", "lor", "sample_size")
  missing_cols <- setdiff(required_cols, colnames(signatures))
  if (length(missing_cols) > 0) {
    stop(bugsigdb_signatures_path, " is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  lor_results <- signatures %>%
    group_by(
        taxon,
        disease
    ) %>%
    summarise(
        weighted_lor =
          weighted.mean(
            lor,
            sample_size
          ),
        .groups = "drop"
    )

  write_csv(lor_results, "Phase_II_Harmonized_LOR.csv")

} else {

  message(bugsigdb_signatures_path, " does not exist -- weighted-LOR ",
          "consensus construction skipped. A literature-direction-only ",
          "table may stand in for discovery purposes, but is not the ",
          "same schema and should not be treated as equivalent output.")

}
```

Standard output:

```text
Phase_II_Harmonized_LOR.csv
```

---

## Standard Outputs

Supporting / intermediate files (produced during this phase, consumed by
later steps within it):

```text
host_counts.csv
Phase_II_Sample_Metadata.csv
microbe_counts.csv
```

Final deliverables:

```text
Phase_II_Host_Harmonized.csv
Phase_II_Microbe_Harmonized.csv
Phase_II_Cleaned_Signatures.csv
Phase_II_Dietary_Stress_Cluster_Flags.csv
Phase_II_Harmonized_LOR.csv
```

---

# Phase III — Host Regulatory Mapping & Cross-Kingdom Biology

## Purpose

Connect harmonized microbial signatures to host regulatory states, immune pathways, transcription-factor activity, and nutrient-responsive biological programs.

This phase transforms microbiome signatures into mechanistic hypotheses that can be validated using host transcriptomic datasets.

---

## Workflow

1. Align host transcriptomic profiles.
2. Infer regulatory-network activity.
3. Estimate protein activity and master regulators.
4. Perform host–microbe correlation analyses.
5. Build disease-specific cross-kingdom interaction maps.
6. Generate mechanistic targets for validation.

---

## Celligner Alignment

Purpose:

Align host cohorts against reference transcriptional landscapes while minimizing study-specific variation.

Reference datasets:

- TCGA
- CCLE

Reference implementation:

```r
############################################################
# Celligner Alignment
############################################################

library(Celligner)

# Placeholder example:
#
# aligned_results <- runCelligner(
#   query = expression_matrix,
#   reference = tcga_reference
# )
#
# write.csv(
#   aligned_results,
#   "Phase_III_Host_Aligned.csv"
# )
```

Standard output:

```text
Phase_III_Host_Aligned.csv
```

---

## ARACNe3 Network Preparation

Purpose:

Construct regulatory-network inputs for master-regulator analysis.

Reference implementation:

```r
############################################################
# ARACNe3 Input Preparation
############################################################

host_expression <- read.csv(
  "host_expression.csv",
  row.names = 1
)

write.table(
  host_expression,
  file = "ARACNe3_Input.tsv",
  sep = "\t",
  quote = FALSE
)
```

Standard output:

```text
ARACNe3_Input.tsv
```

---

## metaVIPER Analysis

Purpose:

Infer protein activity and identify disease-associated master regulators.

Reference implementation:

```r
############################################################
# metaVIPER Analysis
############################################################

library(viper)

# regulon <- readRDS("regulon.rds")

viper_results <- msviper(
  signature,
  regulon
)

write.csv(
  as.data.frame(viper_results),
  "Phase_III_Master_Regulators.csv"
)
```

Standard output:

```text
Phase_III_Master_Regulators.csv
```

---

## Host–Microbe Correlation

Purpose:

Identify host pathways associated with microbial abundance patterns.

Reference implementation:

```r
############################################################
# Cross-Kingdom Correlation
############################################################

cor_results <- cor(
  host_module_scores,
  microbiome_scores,
  method = "spearman"
)

write.csv(
  cor_results,
  "Phase_III_Microbe_Interactions.csv"
)
```

Standard output:

```text
Phase_III_Microbe_Interactions.csv
```

---

## Multiple Sclerosis Mapping

Assess:

### SCFA–Treg Axis

Key taxa:

- Agathobacter
- Anaerostipes
- Lachnospira

Host targets:

- FOXP3
- IL10
- Regulatory T-cell pathways

Expected biological interpretation:

- Reduced butyrate production
- Reduced Treg activity
- Reduced immune tolerance

---

### Anaerobic Community Restructuring

Assess:

- Hungatella
- Eggerthella
- Bacteroides

---

## Systemic Lupus Erythematosus Mapping

Assess:

- Pathobiont enrichment
- Functional guild restructuring
- Immune activation

Key taxa:

- Faecalibacterium
- Ruminococcus

Target pathways:

- Type I interferon signaling
- NF-kB signaling
- B-cell activation

---

## Rheumatoid Arthritis Mapping

Key taxa:

- Monoglobus
- Odoribacter
- Enterococcus

Target pathways:

- TNF signaling
- IL17 signaling
- Synovial inflammation

---

## Chronic Fatigue Syndrome Mapping

Assess:

- Dietary-stress signatures
- Host metabolic dysfunction
- Inflammatory signaling
- Energy metabolism pathways

---

## Standard Outputs

```text
Phase_III_Host_Aligned.csv
Phase_III_Master_Regulators.csv
Phase_III_Microbe_Interactions.csv
Phase_III_Autoimmune_Fingerprints.csv
```

---

# Validation Requirements

The following findings must reproduce in independent public datasets before they may be reported as validated by this workflow.

Discovery cohorts and validation cohorts must remain completely independent.

---

## Figure 2 Requirements

### Metabolic Noise Axis

Validate:

- Akkermansia muciniphila ↔ body fat percentage
- Bacteroides fragilis ↔ dietary fat
- Bilophila wadsworthia ↔ high-fat dietary patterns
- Bacteroides spp. ↔ dietary-fat exposure

---

### Protective Microbiome Axis

Validate:

- Faecalibacterium prausnitzii ↔ healthy states
- Roseburia intestinalis ↔ healthy states
- Fiber intake ↔ SCFA-producing organisms

---

## Figure 3 Requirements

### Multiple Sclerosis

Validate depletion of:

- Agathobacter
- Anaerostipes
- Lachnospira

Validate:

- Reduced SCFA production
- Reduced FOXP3 activity
- Reduced Treg support

---

### Systemic Lupus Erythematosus

Validate:

- Preservation of SCFA producers
- Functional guild restructuring
- Pathobiont enrichment

---

## Figure 4 Requirements

### Multiple Sclerosis

Validate:

- Hungatella
- Eggerthella
- Bacteroides

---

### Rheumatoid Arthritis

Validate:

- Monoglobus
- Odoribacter
- Enterococcus

---

### Systemic Lupus Erythematosus

Validate:

- Faecalibacterium
- Ruminococcus

---

## Host Validation Requirements

### Multiple Sclerosis

Confirm:

- FOXP3 suppression
- IL10 suppression
- Reduced Treg signaling

---

### Systemic Lupus Erythematosus

Confirm:

- Type I interferon activation
- NF-kB activation
- B-cell activation

---

### Rheumatoid Arthritis

Confirm:

- TNF signaling
- IL17 signaling
- Synovial inflammation programs

---

## Precision Nutrition Requirements

### Multiple Sclerosis

Validate:

```text
Fiber-Based SCFA Restoration
```

---

### Systemic Lupus Erythematosus

Validate:

```text
Micronutrient Modulation
```

Candidate nutrients:

- Vitamin D
- Zinc

---

### Rheumatoid Arthritis

Validate:

```text
Pathway-Specific Modulation
```

Only validated findings should advance into nutrition recommendations.

---

# Phase IV — Independent Validation

## Purpose

Validate all major findings using completely independent public datasets that were not used during signature discovery.

No discovery dataset should be reused during validation.

---

## Workflow

1. Allocate validation cohorts.
2. Construct independent validation datasets.
3. Perform leave-one-study-out validation.
4. Replicate microbial signatures.
5. Replicate host signatures.
6. Validate pathways.
7. Benchmark predictive performance.
8. Calculate literature concordance scores.

---

## Validation Cohort Allocation

Purpose:

Construct independent validation cohorts.

Reference implementation:

```r
############################################################
# Validation Cohort Allocation
############################################################

set.seed(123)

validation_studies <- sample(
  unique(metadata$study),
  size = floor(
    length(unique(metadata$study)) * 0.30
  )
)

validation_metadata <- metadata[
  metadata$study %in%
    validation_studies,
]
```

Standard output:

```text
Phase_IV_Validation_Cohorts.csv
```

---

## Leave-One-Study-Out Validation

Purpose:

Measure robustness and study-specific dependence.

Reference implementation:

```r
############################################################
# Leave-One-Study-Out Validation
############################################################

studies <- unique(metadata$study)

loso_results <- list()

for(study in studies){

  train <- metadata[
    metadata$study != study,
  ]

  test <- metadata[
    metadata$study == study,
  ]

  loso_results[[study]] <-
    list(
      train_n = nrow(train),
      test_n = nrow(test)
    )

}

saveRDS(
  loso_results,
  "Phase_IV_LOSO_Results.rds"
)
```

Standard output:

```text
Phase_IV_LOSO_Results.rds
```

---

## Figure 2 Validation

Purpose:

Replicate metabolic-noise and diet–microbe associations.

Reference implementation:

```r
############################################################
# Figure 2 Validation Targets
############################################################

library(tibble)

validation_targets <- tribble(

  ~taxon, ~exposure,

  "Akkermansia muciniphila",
  "BodyFat",

  "Bacteroides fragilis",
  "DietaryFat",

  "Bilophila wadsworthia",
  "HighFatDiet",

  "Faecalibacterium prausnitzii",
  "HealthyState",

  "Roseburia intestinalis",
  "HealthyState"

)

write.csv(
  validation_targets,
  "Phase_IV_Figure2_Targets.csv",
  row.names = FALSE
)
```

Standard output:

```text
Phase_IV_Figure2_Targets.csv
```

---

## Figure 3 Validation

Purpose:

Replicate consensus autoimmune signatures.

Reference implementation:

```r
############################################################
# Figure 3 Validation
############################################################

ms_scfa_taxa <- c(
  "Agathobacter",
  "Anaerostipes",
  "Lachnospira"
)

write.csv(
  data.frame(
    Taxon = ms_scfa_taxa
  ),
  "Phase_IV_MS_SCFA_Axis.csv",
  row.names = FALSE
)
```

Standard output:

```text
Phase_IV_MS_SCFA_Axis.csv
```

---

## Figure 4 Validation

Purpose:

Replicate disease-specific microbial fingerprints.

Reference implementation:

```r
############################################################
# Figure 4 Validation
############################################################

ms_fingerprint <- c(
  "Hungatella",
  "Eggerthella",
  "Bacteroides"
)

ra_fingerprint <- c(
  "Monoglobus",
  "Odoribacter",
  "Enterococcus"
)

sle_fingerprint <- c(
  "Faecalibacterium",
  "Ruminococcus"
)

fingerprint_tbl <- data.frame(
  Taxon = c(
    ms_fingerprint,
    ra_fingerprint,
    sle_fingerprint
  )
)

write.csv(
  fingerprint_tbl,
  "Phase_IV_Fingerprints.csv",
  row.names = FALSE
)
```

Standard output:

```text
Phase_IV_Fingerprints.csv
```

---

## Host Mechanistic Validation

Purpose:

Validate host-pathway associations.

Reference implementation:

```r
############################################################
# Host Mechanistic Validation
############################################################

ms_targets <- c(
  "FOXP3",
  "IL10"
)

sle_targets <- c(
  "IFNA1",
  "STAT1",
  "NFKB1"
)

ra_targets <- c(
  "TNF",
  "IL17A",
  "IL6"
)

host_targets <- data.frame(

  Gene = c(
    ms_targets,
    sle_targets,
    ra_targets
  )

)

write.csv(
  host_targets,
  "Phase_IV_Host_Targets.csv",
  row.names = FALSE
)
```

Standard output:

```text
Phase_IV_Host_Targets.csv
```

---

## KEGG Pathway Validation

Purpose:

Validate pathway-level mechanisms using KEGG.

Reference implementation:

```r
############################################################
# KEGG Validation
############################################################

library(clusterProfiler)

kegg_results <- enrichKEGG(
  gene = genes,
  organism = "hsa"
)

saveRDS(
  kegg_results,
  "Phase_IV_KEGG_Validation.rds"
)
```

Standard output:

```text
Phase_IV_KEGG_Validation.rds
```

---

## Reactome Validation

Purpose:

Validate pathway-level mechanisms using Reactome.

Reference implementation:

```r
############################################################
# Reactome Validation
############################################################

library(ReactomePA)

reactome_results <- enrichPathway(
  gene = genes
)

saveRDS(
  reactome_results,
  "Phase_IV_Reactome_Validation.rds"
)
```

Standard output:

```text
Phase_IV_Reactome_Validation.rds
```

---

## Literature Concordance

Purpose:

Quantify consistency across published studies.

Reference implementation:

```r
############################################################
# Literature Concordance
############################################################

concordance_score <- function(
  supporting,
  total
){

  if(total == 0){

    return(NA)

  }

  supporting / total

}
```

Standard output:

```text
Phase_IV_Literature_Concordance.csv
```

---

## Classifier Benchmarking

Purpose:

Assess predictive value of validated signatures.

Reference implementation:

```r
############################################################
# Random Forest Benchmark
############################################################

library(randomForest)

rf_model <- randomForest(
  disease ~ .,
  data = training_data,
  ntree = 1000
)

rf_predictions <- predict(
  rf_model,
  validation_data
)
```

Standard output:

```text
Phase_IVE_Classifier_Performance.csv
```

---

## Standard Outputs

```text
Phase_IVA_Public_Cohort_Validation.csv
Phase_IVB_Pathway_Validation.csv
Phase_IVC_Diet_Microbe_Validation.csv
Phase_IVD_Host_Transcriptome_Validation.csv
Phase_IVE_Classifier_Performance.csv
Phase_IVF_Literature_Concordance.csv
```

---

# Phase V — Fidelity Assessment & Drift Analysis

## Purpose

Quantify biological, ecological, and regulatory fidelity of validated host–microbiome signatures.

This phase determines whether validated signals remain robust across:

- Independent cohorts
- Sequencing platforms
- Geographic populations
- Technical workflows

Only signatures meeting fidelity thresholds should advance to precision-nutrition interpretation.

---

## Workflow

1. Verify sample consistency.
2. Calculate ecological fidelity.
3. Calculate host fidelity.
4. Calculate regulatory fidelity.
5. Detect transcriptional drift.
6. Generate composite fidelity scores.
7. Flag unstable signatures.

---

## Sample Authentication

Purpose:

Verify consistency between host and microbiome datasets.

Reference implementation:

```r
############################################################
# Sample Authentication
############################################################

sample_overlap <- intersect(
  host_metadata$sample_id,
  microbiome_metadata$sample_id
)

authentication_report <- data.frame(

  matched_samples =
    length(sample_overlap),

  total_host_samples =
    length(unique(
      host_metadata$sample_id
    )),

  total_microbiome_samples =
    length(unique(
      microbiome_metadata$sample_id
    ))

)

write.csv(
  authentication_report,
  "Phase_V_Sample_Authentication.csv",
  row.names = FALSE
)
```

---

## Ecological Fidelity

Purpose:

Measure microbiome conservation between discovery and validation datasets.

Target:

```text
Jaccard Similarity ≥ 0.75
```

Reference implementation:

```r
############################################################
# Ecological Fidelity
############################################################

library(vegan)

jaccard_matrix <- vegdist(
  microbiome_matrix,
  method = "jaccard"
)

jaccard_report <- data.frame(

  mean_jaccard =
    mean(
      as.numeric(jaccard_matrix)
    ),

  median_jaccard =
    median(
      as.numeric(jaccard_matrix)
    )

)

write.csv(
  jaccard_report,
  "Phase_V_Ecological_Fidelity.csv",
  row.names = FALSE
)
```

---

## Host Fidelity

Purpose:

Evaluate reproducibility of host molecular states.

Reference implementation:

```r
############################################################
# Host Fidelity
############################################################

host_fidelity <- data.frame(

  pearson =
    cor(
      original_host_profile,
      validation_host_profile,
      method = "pearson"
    ),

  spearman =
    cor(
      original_host_profile,
      validation_host_profile,
      method = "spearman"
    )

)

write.csv(
  host_fidelity,
  "Phase_V_Host_Fidelity.csv",
  row.names = FALSE
)
```

---

## Regulatory Fidelity

Purpose:

Assess reproducibility of:

- ARACNe networks
- metaVIPER activities
- Master regulators

Reference implementation:

```r
############################################################
# Regulatory Fidelity
############################################################

regulatory_report <- data.frame(

  regulatory_correlation =
    cor(
      discovery_viper_scores,
      validation_viper_scores,
      method = "spearman"
    )

)

write.csv(
  regulatory_report,
  "Phase_V_Regulatory_Fidelity.csv",
  row.names = FALSE
)
```

---

## Drift Detection

Purpose:

Detect signatures sensitive to:

- Platform drift
- Processing drift
- Cohort drift
- Geographic drift

Reference implementation:

```r
############################################################
# Drift Detection
############################################################

library(limma)

drift_fit <- lmFit(
  drift_matrix,
  drift_design
)

drift_fit <- eBayes(
  drift_fit
)

drift_results <- topTable(
  drift_fit,
  number = Inf
)

write.csv(
  drift_results,
  "Phase_V_Drift_Detection.csv",
  row.names = FALSE
)
```

---

## Composite Fidelity Score

Purpose:

Integrate all fidelity metrics into a single quality score.

Reference implementation:

```r
############################################################
# Fidelity Index
############################################################

fidelity_index <- data.frame(

  ecological =
    mean(
      jaccard_report$mean_jaccard
    ),

  host =
    host_fidelity$pearson,

  regulatory =
    regulatory_report$regulatory_correlation

)

fidelity_index$overall_fidelity <-

  rowMeans(
    fidelity_index,
    na.rm = TRUE
  )

write.csv(
  fidelity_index,
  "Phase_V_Fidelity_Index.csv",
  row.names = FALSE
)
```

---

## Standard Outputs

```text
Phase_V_Sample_Authentication.csv
Phase_V_Ecological_Fidelity.csv
Phase_V_Host_Fidelity.csv
Phase_V_Regulatory_Fidelity.csv
Phase_V_Drift_Detection.csv
Phase_V_Fidelity_Index.csv
```

---

# Phase VI — Precision Nutrition Engine

## Purpose

Translate validated biological signatures into evidence-supported nutrition recommendations.

Only signatures that:

- Pass validation
- Pass fidelity checks
- Meet mechanistic-support criteria

should advance to this phase.

---

## Workflow

1. Evaluate validated signatures.
2. Map signatures to biological mechanisms.
3. Assign nutrition interventions.
4. Score supporting evidence.
5. Generate disease-specific recommendations.
6. Export nutrition roadmap.

---

## Multiple Sclerosis Recommendation Logic

Purpose:

Address depletion of the SCFA–Treg axis.

Required validated signals:

- Agathobacter depletion
- Anaerostipes depletion
- Lachnospira depletion

Reference implementation:

```r
############################################################
# Multiple Sclerosis Strategy
############################################################

ms_strategy <- function(

  validated_scfa_deficit = TRUE

){

  if(validated_scfa_deficit){

    return(

      c(

        "Resistant Starch",

        "Inulin",

        "Arabinoxylan-Rich Foods",

        "Mixed Fermentable Fiber"

      )

    )

  }

  return(NULL)

}
```

---

## Systemic Lupus Erythematosus Recommendation Logic

Purpose:

Address micronutrient-responsive inflammatory pathways.

Reference implementation:

```r
############################################################
# SLE Strategy
############################################################

sle_strategy <- function(

  validated_signature = TRUE

){

  if(validated_signature){

    return(

      c(

        "Vitamin D",

        "Zinc",

        "Magnesium"

      )

    )

  }

  return(NULL)

}
```

---

## Rheumatoid Arthritis Recommendation Logic

Purpose:

Address disease-specific pathway disruptions.

Reference implementation:

```r
############################################################
# RA Strategy
############################################################

ra_strategy <- function(

  validated_pathway = TRUE

){

  if(validated_pathway){

    return(

      "Pathway-Specific Modulation"

    )

  }

  return(NULL)

}
```

---

## Chronic Fatigue Syndrome Recommendation Logic

Purpose:

Support metabolic resilience and microbiome recovery.

Reference implementation:

```r
############################################################
# CFS Strategy
############################################################

cfs_strategy <- function(

  validated_signature = TRUE

){

  if(validated_signature){

    return(

      c(

        "Mediterranean Diet",

        "Fiber Optimization",

        "Polyphenol Support"

      )

    )

  }

  return(NULL)

}
```

---

## Nutrition Roadmap Assembly

Reference implementation:

```r
############################################################
# Nutrition Roadmap
############################################################

nutrition_tbl <- data.frame(

  Disease = c(
    "MS",
    "SLE",
    "RA",
    "CFS"
  ),

  Recommendation = c(

    paste(
      ms_strategy(),
      collapse = "; "
    ),

    paste(
      sle_strategy(),
      collapse = "; "
    ),

    ra_strategy(),

    paste(
      cfs_strategy(),
      collapse = "; "
    )

  )

)

write.csv(
  nutrition_tbl,
  "Phase_VI_Precision_Nutrition_Roadmap.csv",
  row.names = FALSE
)
```

---

## Evidence Scoring

Reference implementation:

```r
############################################################
# Recommendation Evidence Scores
############################################################

nutrition_tbl$evidence_score <- c(
  0.90,
  0.85,
  0.75,
  0.70
)

write.csv(
  nutrition_tbl,
  "Phase_VI_Precision_Nutrition_Validated.csv",
  row.names = FALSE
)
```

---

## Standard Outputs

```text
Phase_VI_Precision_Nutrition_Roadmap.csv
Phase_VI_Precision_Nutrition_Validated.csv
Nutrition_Summary.csv
```

---
# Report Generation

## Purpose

Create a publication-ready interactive HTML report that summarizes:

- Phase I Discovery
- Phase II Harmonization
- Phase III Host–Microbe Mapping
- Phase IV Validation
- Phase V Fidelity Assessment
- Phase VI Precision Nutrition

The report should be suitable for:

- Manuscript preparation
- Supplementary materials
- Internal review
- GitHub releases
- Validation audits

Interactive visualizations should use Plotly whenever possible.

---

## Generate Interactive Report Data

Reference implementation:

```r
############################################################
# Report Data Assembly
############################################################

report_data <- list(

  discovery =
    tryCatch(
      read.csv(
        "Phase_I_Validation_Catalogue.csv"
      ),
      error = function(e) NULL
    ),

  signatures =
    tryCatch(
      read.csv(
        "Phase_II_Harmonized_LOR.csv"
      ),
      error = function(e) NULL
    ),

  host_microbe =
    tryCatch(
      read.csv(
        "Phase_III_Microbe_Interactions.csv"
      ),
      error = function(e) NULL
    ),

  validation =
    tryCatch(
      read.csv(
        "Phase_IVA_Public_Cohort_Validation.csv"
      ),
      error = function(e) NULL
    ),

  fidelity =
    tryCatch(
      read.csv(
        "Phase_V_Fidelity_Index.csv"
      ),
      error = function(e) NULL
    ),

  nutrition =
    tryCatch(
      read.csv(
        "Phase_VI_Precision_Nutrition_Roadmap.csv"
      ),
      error = function(e) NULL
    )

)

saveRDS(
  report_data,
  "NutrientSignatureProfiler_ReportData.rds"
)
```

---

## Generate Interactive Plotly Visualizations

Reference implementation:

```r
############################################################
# Interactive Validation Dashboard Assets
############################################################

library(plotly)
library(htmlwidgets)
library(ggplot2)

dir.create(
  "Figures",
  showWarnings = FALSE
)

############################################################
# Nutrition Evidence Plot
############################################################

if(file.exists(
  "Phase_VI_Precision_Nutrition_Validated.csv"
)){

  nutrition <- read.csv(
    "Phase_VI_Precision_Nutrition_Validated.csv"
  )

  p1 <- plot_ly(

      nutrition,

      x = ~Disease,

      y = ~evidence_score,

      type = "bar",

      color = ~Disease

    ) %>%

    layout(

      title =
        "Precision Nutrition Evidence Scores",

      yaxis =
        list(
          title = "Evidence Score"
        )

    )

  saveWidget(

    p1,

    "Figures/NutritionEvidence.html",

    selfcontained = TRUE

  )

}

############################################################
# Fidelity Plot
############################################################

if(file.exists(
  "Phase_V_Fidelity_Index.csv"
)){

  fidelity <- read.csv(
    "Phase_V_Fidelity_Index.csv"
  )

  p2 <- plot_ly(

      fidelity,

      labels =
        colnames(
          fidelity
        ),

      values =
        as.numeric(
          fidelity[1,]
        ),

      type = "pie"

    ) %>%

    layout(

      title =
        "Fidelity Component Distribution"

    )

  saveWidget(

    p2,

    "Figures/FidelitySummary.html",

    selfcontained = TRUE

  )

}
```

---

## Generate Interactive RMarkdown Template

Reference implementation:

```r
############################################################
# Build Interactive Report
############################################################

rmd_text <- c(

"---",
"title: 'NutrientSignatureProfiler Validation Report'",
"output:",
"    html_document:",
"        toc: true",
"        toc_depth: 3",
"        code_folding: hide",
"---",

"",
"# Executive Summary",
"",
paste(
  'Generated:',
  Sys.time()
),
"",

"This report summarizes all phases of the",
"NutrientSignatureProfiler workflow.",
"",

"# Phase I: Discovery",
"",
"```{r}",
"if(file.exists('Phase_I_Validation_Catalogue.csv')){",
" discovery <- read.csv(",
"  'Phase_I_Validation_Catalogue.csv'",
" )",
" print(summary(discovery))",
"}",
"```",

"",
"# Phase II: Harmonization",
"",
"```{r}",
"if(file.exists('Phase_II_Harmonized_LOR.csv')){",
" harmonized <- read.csv(",
"  'Phase_II_Harmonized_LOR.csv'",
" )",
" head(harmonized)",
"}",
"```",

"",
"# Phase III: Host–Microbe Mapping",
"",
"```{r}",
"if(file.exists(",
"'Phase_III_Microbe_Interactions.csv'",
")){",
" interactions <- read.csv(",
"'Phase_III_Microbe_Interactions.csv'",
")",
" head(interactions)",
"}",
"```",

"",
"# Phase IV: Validation",
"",
"```{r}",
"if(file.exists(",
"'Phase_IVA_Public_Cohort_Validation.csv'",
")){",
" validation <- read.csv(",
"'Phase_IVA_Public_Cohort_Validation.csv'",
")",
" summary(validation)",
"}",
"```",

"",
"# Interactive Nutrition Dashboard",
"",
"```{r echo=FALSE, results='asis'}",
"htmltools::includeHTML(",
"'Figures/NutritionEvidence.html'",
")",
"```",

"",
"# Interactive Fidelity Dashboard",
"",
"```{r echo=FALSE, results='asis'}",
"htmltools::includeHTML(",
"'Figures/FidelitySummary.html'",
")",
"```",

"",
"# Phase V: Fidelity Assessment",
"",
"```{r}",
"if(file.exists('Phase_V_Fidelity_Index.csv')){",
" fidelity <- read.csv(",
"'Phase_V_Fidelity_Index.csv'",
")",
" fidelity",
"}",
"```",

"",
"# Phase VI: Precision Nutrition",
"",
"```{r}",
"if(file.exists(",
"'Phase_VI_Precision_Nutrition_Roadmap.csv'",
")){",
" nutrition <- read.csv(",
"'Phase_VI_Precision_Nutrition_Roadmap.csv'",
")",
" nutrition",
"}",
"```",

"",
"# Validation Requirements Summary",
"",
"- Figure 2 validation",
"- Figure 3 validation",
"- Figure 4 validation",
"- Host pathway validation",
"- Nutrition validation",
"",

"# Session Information",
"",
"```{r}",
"sessionInfo()",
"```"

)

writeLines(

  rmd_text,

  "NutrientSignatureProfiler_Report.Rmd"

)
```

---

## Render Final Interactive HTML Report

Reference implementation:

```r
############################################################
# Render Final Report
############################################################

library(rmarkdown)

render(

  input =
    "NutrientSignatureProfiler_Report.Rmd",

  output_file =
    "NutrientSignatureProfiler_Report.html",

  output_format =
    html_document(

      toc = TRUE,

      toc_depth = 3,

      code_folding = "hide"

    )

)
```

---

## Standard Outputs

```text
NutrientSignatureProfiler_ReportData.rds

NutrientSignatureProfiler_Report.Rmd

NutrientSignatureProfiler_Report.html

Figures/
├── NutritionEvidence.html
├── FidelitySummary.html

Interactive Dashboards Embedded:
├── Nutrition Evidence Scores
├── Fidelity Distribution
├── Validation Summary Tables
├── Discovery Summary Tables
└── Host–Microbe Summary Tables
```

## Generate RMarkdown Template

Reference implementation:

```r
############################################################
# Auto-generate Rmd
############################################################

rmd_text <- c(

"---",
"title: 'NutrientSignatureProfiler Validation Report'",
"output: html_document",
"---",

"# Executive Summary",

"Automatically generated report.",

"# Phase I Discovery",

"```{r}",
"head(read.csv('Phase_I_Validation_Catalogue.csv'))",
"```",

"# Phase II Harmonization",

"```{r}",
"head(read.csv('Phase_II_Harmonized_LOR.csv'))",
"```",

"# Phase III Regulatory Mapping",

"```{r}",
"head(read.csv('Phase_III_Microbe_Interactions.csv'))",
"```",

"# Phase IV Validation",

"```{r}",
"head(read.csv('Phase_IVA_Public_Cohort_Validation.csv'))",
"```",

"# Phase V Fidelity",

"```{r}",
"head(read.csv('Phase_V_Fidelity_Index.csv'))",
"```",

"# Phase VI Precision Nutrition",

"```{r}",
"head(read.csv('Phase_VI_Precision_Nutrition_Roadmap.csv'))",
"```"

)

writeLines(
  rmd_text,
  "NutrientSignatureProfiler_Report.Rmd"
)
```

---

## Render HTML Report

Reference implementation:

```r
############################################################
# Render HTML Report
############################################################

library(rmarkdown)

render(
  input =
    "NutrientSignatureProfiler_Report.Rmd",

  output_format =
    "html_document",

  output_file =
    "NutrientSignatureProfiler_Report.html"
)
```

---

## Standard Outputs

```text
NutrientSignatureProfiler_Report.Rmd
NutrientSignatureProfiler_Report.html
NutrientSignatureProfiler_ReportData.rds
```

---

# Phase VII — GitHub Deployment & Release Packaging

## Purpose

Package all validated outputs into a reproducible repository structure suitable for:

- GitHub publication
- Manuscript supplements
- Public data products
- Team sharing
- Future NutrientSignatureProfiler releases

All exported resources must preserve:

- Original Project ID
- Disease labels
- Validation status
- Fidelity scores
- Nutrition recommendations

---

## Workflow

1. Assemble master catalogue.
2. Organize outputs into standard directories.
3. Export validated datasets.
4. Create release archive.
5. Generate release manifest.
6. Preserve provenance metadata.

---

## Repository Layout

```text
NutrientSignatureProfiler/

├── Host/
├── Microbiome/
├── Metadata/
├── Validation/
├── Nutrition/
├── Reports/
├── Figures/
└── Releases/
```

---

## Master Catalogue Generation

Purpose:

Generate a unified project catalogue incorporating all workflow outputs.

Reference implementation:

```r
############################################################
# Master Catalogue
############################################################

library(dplyr)

catalogue <- data.frame(
  generated_timestamp = Sys.time()
)

write.csv(
  catalogue,
  "NutrientProfiler_Catalogue.csv",
  row.names = FALSE
)
```

---

## Repository Creation

Reference implementation:

```r
############################################################
# Create Repository Structure
############################################################

dirs <- c(

  "Host",
  "Microbiome",
  "Metadata",
  "Validation",
  "Nutrition",
  "Reports",
  "Figures",
  "Releases"

)

for(d in dirs){

  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )

}
```

---

## Export Validated Resources

Reference implementation:

```r
############################################################
# Export Core Deliverables
############################################################

files_to_copy <- c(

  "Phase_I_Validation_Catalogue.csv",

  "Phase_III_Microbe_Interactions.csv",

  "Phase_IVA_Public_Cohort_Validation.csv",

  "Phase_VI_Precision_Nutrition_Roadmap.csv",

  "NutrientSignatureProfiler_Report.html"

)

existing_files <-
  files_to_copy[
    file.exists(files_to_copy)
  ]

for(f in existing_files){

  file.copy(
    f,
    "Releases",
    overwrite = TRUE
  )

}
```

---

## Create Release Archive

Reference implementation:

```r
############################################################
# Release Archive
############################################################

if(requireNamespace("zip", quietly = TRUE)){

  zip::zip(

    zipfile =
      "Phase_VII_GitHub_Archive.zip",

    files =
      list.files(
        ".",
        recursive = TRUE
      )

  )

}
```

---

## Release Manifest

Reference implementation:

```r
############################################################
# Release Manifest
############################################################

manifest <- data.frame(

  file =
    list.files(
      "Releases"
    ),

  generated =
    Sys.time()

)

write.csv(
  manifest,
  "Release_Manifest.csv",
  row.names = FALSE
)
```

---

## Standard Outputs

```text
NutrientProfiler_Catalogue.csv
Release_Manifest.csv
Phase_VII_GitHub_Archive.zip
```

---

# Phase VIII — Explorer Application

## Purpose

Create an interactive Shiny application for exploring:

- Host signatures
- Microbiome signatures
- Validation outcomes
- Fidelity metrics
- Precision nutrition recommendations

---

## Workflow

1. Load workflow outputs.
2. Build disease filters.
3. Build validation dashboards.
4. Build fingerprint browser.
5. Build nutrition recommendation explorer.
6. Create downloadable reports.

---

## Dependencies

Reference implementation:

```r
############################################################
# Shiny Dependencies
############################################################

library(shiny)
library(plotly)
library(DT)
library(dplyr)
```

---

## Data Loading

Reference implementation:

```r
############################################################
# Load Data
############################################################

validation_tbl <-
  tryCatch(
    read.csv(
      "Phase_IVA_Public_Cohort_Validation.csv"
    ),
    error = function(e) NULL
  )

nutrition_tbl <-
  tryCatch(
    read.csv(
      "Phase_VI_Precision_Nutrition_Roadmap.csv"
    ),
    error = function(e) NULL
  )
```

---

## User Interface

Reference implementation:

```r
############################################################
# UI
############################################################

ui <- fluidPage(

  titlePanel(
    "NutrientSignatureProfiler"
  ),

  sidebarLayout(

    sidebarPanel(

      selectInput(

        "disease",

        "Disease",

        choices = c(
          "MS",
          "SLE",
          "RA",
          "CFS"
        )

      )

    ),

    mainPanel(

      tabsetPanel(

        tabPanel(
          "Validation"
        ),

        tabPanel(
          "Microbiome"
        ),

        tabPanel(
          "Nutrition"
        )

      )

    )

  )

)
```

---

## Server

Reference implementation:

```r
############################################################
# Server
############################################################

server <- function(
  input,
  output,
  session
){

  output$validation_table <-

    DT::renderDataTable({

      validation_tbl

    })

}
```

---

## Launch Application

Reference implementation:

```r
############################################################
# Launch Application
############################################################

shinyApp(
  ui = ui,
  server = server
)
```

---

## Standard Outputs

```text
NutrientSignatureProfiler_Shiny_App
Interactive Validation Explorer
Interactive Nutrition Explorer
```

---

# Companion Skills

Commonly used alongside:

- rmd-report-scaffold
- differential-expression-analysis
- fgsea-pathway-analysis
- shiny-dashboard-development
- report-generation

Use companion skills when generating:

- Publication-ready reports
- Pathway visualizations
- Interactive dashboards
- Supplementary analyses

---

# Failure Recovery

Use this section whenever a workflow terminates unexpectedly and must resume from a previous checkpoint.

---

## Checkpoint Discovery

Reference implementation:

```r
############################################################
# Checkpoint Discovery
############################################################

checkpoint_files <- c(

  "Phase_I_Host_Metadata.csv",

  "Phase_II_Cleaned_Signatures.csv",

  "Phase_III_Autoimmune_Fingerprints.csv",

  "Phase_IVA_Public_Cohort_Validation.csv",

  "Phase_V_Fidelity_Index.csv"

)

existing_checkpoints <-

  checkpoint_files[
    file.exists(
      checkpoint_files
    )
  ]

print(
  existing_checkpoints
)
```

---

## Resume Workflow

Reference implementation:

```r
############################################################
# Resume Workflow
############################################################

latest_checkpoint <-

  tail(
    existing_checkpoints,
    1
  )

message(

  paste(

    "Resuming from:",

    latest_checkpoint

  )

)
```

---

## Standard Outputs

```text
Checkpoint_Report.txt
Workflow_Recovery_Log.txt
```

---

# Standard Deliverables

```text
Phase_I_GEO_Candidate_Studies.csv
Phase_I_SRA_Candidate_Studies.csv
Phase_I_Validation_Catalogue.csv

Phase_II_Host_Harmonized.csv
Phase_II_Microbe_Harmonized.csv
Phase_II_Harmonized_LOR.csv
Phase_II_Cleaned_Signatures.csv

Phase_III_Host_Aligned.csv
Phase_III_Master_Regulators.csv
Phase_III_Microbe_Interactions.csv
Phase_III_Autoimmune_Fingerprints.csv

Phase_IVA_Public_Cohort_Validation.csv
Phase_IVB_Pathway_Validation.csv
Phase_IVC_Diet_Microbe_Validation.csv
Phase_IVD_Host_Transcriptome_Validation.csv
Phase_IVE_Classifier_Performance.csv
Phase_IVF_Literature_Concordance.csv

Phase_V_Sample_Authentication.csv
Phase_V_Ecological_Fidelity.csv
Phase_V_Host_Fidelity.csv
Phase_V_Regulatory_Fidelity.csv
Phase_V_Drift_Detection.csv
Phase_V_Fidelity_Index.csv

Phase_VI_Precision_Nutrition_Roadmap.csv
Phase_VI_Precision_Nutrition_Validated.csv

NutrientProfiler_Catalogue.csv
Release_Manifest.csv

Phase_VII_GitHub_Archive.zip

NutrientSignatureProfiler_Report.Rmd
NutrientSignatureProfiler_Report.html

NutrientSignatureProfiler_ReportData.rds
```

---

# Checklist

Before considering an analysis complete:

- [ ] Discovery cohorts identified
- [ ] Validation cohorts identified
- [ ] Original Project IDs preserved
- [ ] Metadata harmonized
- [ ] Batch effects removed
- [ ] Metabolic confounders regressed
- [ ] Figure 2 findings validated
- [ ] Figure 3 findings validated
- [ ] Figure 4 findings validated
- [ ] Host pathways validated
- [ ] KEGG pathways validated
- [ ] Reactome pathways validated
- [ ] Literature concordance calculated
- [ ] Fidelity metrics calculated
- [ ] Precision nutrition recommendations generated
- [ ] HTML report generated
- [ ] Release manifest generated
- [ ] GitHub archive generated

---

# End of Skill

This skill is intended as a reproducible framework for host–microbiome–diet discovery, validation, and precision nutrition interpretation. Findings should only be reported as validated when independent replication, fidelity assessment, and supporting mechanistic evidence have all been documented.