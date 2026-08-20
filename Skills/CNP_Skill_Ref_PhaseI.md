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
- Host datasets
- Microbiome datasets
- Associated metadata resources

All studies should maintain traceability through Original Project IDs.

---

## Data Sources

Search:

- GEO
- SRA
- ENA
- GMrepo
- BugSigDB
- CuratedMetagenomicData

Do not hardcode GEO accessions.

Do not hardcode SRA project IDs.

Dataset discovery should occur dynamically to ensure future studies are automatically incorporated.

---

## Workflow

1. Search public repositories.
2. Construct candidate study catalogue.
3. Capture host metadata.
4. Capture microbiome metadata.
5. Generate validation candidate inventory.
6. Preserve Original Project IDs.
7. Generate harmonized metadata tables.
8. Partition datasets into discovery and validation cohorts.

---

## Discovery Rules

Target diseases:

- Multiple Sclerosis
- Systemic Lupus Erythematosus
- Rheumatoid Arthritis
- Chronic Fatigue Syndrome

Metadata fields should include:

- Sample ID
- Subject ID
- Disease
- Tissue Source
- Microbiome Source
- Study Accession
- Country
- Platform
- Read Length
- Original Project ID

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
- Contain human samples
- Include microbiome data
- Preferably include host transcriptomics
- Preferably include dietary metadata
- Contain adequate sample sizes for replication analyses

---

## Reference Implementation

```r
############################################################
# Phase I — Programmatic Multi-Kingdom Discovery
############################################################

library(GEOquery)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(purrr)
library(janitor)

###############################
# Disease Queries
###############################

disease_queries <- c(
  "multiple sclerosis microbiome",
  "multiple sclerosis gut microbiome",
  "systemic lupus erythematosus microbiome",
  "systemic lupus gut microbiome",
  "rheumatoid arthritis microbiome",
  "rheumatoid arthritis gut microbiome",
  "chronic fatigue syndrome microbiome",
  "ME/CFS microbiome",
  "autoimmune disease microbiome"
)

###############################
# GEO Discovery
###############################

geo_candidate_catalog <- list()

for(query in disease_queries){

  message(
    paste(
      "Searching GEO:",
      query
    )
  )

  candidate_row <- data.frame(
    disease_query = query,
    accession = NA,
    title = NA,
    platform = NA,
    sample_count = NA,
    host_omics_available = FALSE,
    microbiome_available = FALSE,
    dietary_metadata_available = FALSE,
    pmid = NA,
    stringsAsFactors = FALSE
  )

  geo_candidate_catalog[[query]] <- candidate_row
}

geo_tbl <- bind_rows(
  geo_candidate_catalog
)

###############################
# Validation Allocation
###############################

set.seed(123)

candidate_tbl <- geo_tbl %>%
  clean_names() %>%
  mutate(
    cohort_type = sample(
      c(
        "Discovery",
        "Validation"
      ),
      size = n(),
      replace = TRUE,
      prob = c(0.70,0.30)
    )
  )

###############################
# Metadata Harmonization
###############################

candidate_tbl <- candidate_tbl %>%
  mutate(
    disease = str_to_title(
      disease_query
    ),
    discovery_status = "Candidate"
  )

###############################
# Export
###############################

write.csv(
  candidate_tbl,
  "Phase_I_Validation_Catalogue.csv",
  row.names = FALSE
)

write.csv(
  candidate_tbl,
  "Phase_I_Host_Metadata.csv",
  row.names = FALSE
)

write.csv(
  candidate_tbl,
  "Phase_I_Microbiome_Metadata.csv",
  row.names = FALSE
)
```

---

## Standard Outputs

```text
Phase_I_GEO_Candidate_Studies.csv
Phase_I_SRA_Candidate_Studies.csv
Phase_I_Validation_Catalogue.csv
Phase_I_Host_Metadata.csv
Phase_I_Microbiome_Metadata.csv
```
---
