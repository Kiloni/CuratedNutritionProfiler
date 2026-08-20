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
