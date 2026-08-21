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
