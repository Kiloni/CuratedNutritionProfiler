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

1. Normalize host expression.
2. Correct batch effects.
3. Harmonize microbiome signatures.
4. Build weighted meta-signatures.
5. Remove metabolic confounders.
6. Construct cleaned signature objects.
7. Export harmonized matrices.

---

## Host TPM Normalization

Purpose:

Generate normalized host transcriptomic matrices suitable for integrative analyses.

Reference implementation:

```r
library(edgeR)

counts <- read.csv(
  "host_counts.csv",
  row.names = 1
)

dge <- DGEList(
  counts
)

tpm_matrix <- cpm(
  dge,
  normalized.lib.sizes = TRUE
)
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

expr_corrected <- removeBatchEffect(
  expression_matrix,
  batch = metadata$study_id
)
```

Standard output:

```text
Phase_II_Host_Harmonized.csv
```

---

## Metabolic Noise Regression

Purpose:

Remove microbiome signals associated with BMI, body fat, and dietary fat before biomarker construction.

Reference implementation:

```r
microbe_adjusted <- list()

for(feature in colnames(microbe_matrix)){

  fit <- lm(
    microbe_matrix[,feature] ~
      metadata$bmi +
      metadata$body_fat +
      metadata$dietary_fat +
      metadata$age +
      metadata$sex
  )

  microbe_adjusted[[feature]] <- resid(fit)

}
```

Standard output:

```text
Phase_II_Cleaned_Signatures.csv
```

---

## BugSigDB Consensus Signature Construction

Purpose:

Construct weighted disease meta-signatures.

Reference implementation:

```r
library(dplyr)

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
        )
  )
```

Standard output:

```text
Phase_II_Harmonized_LOR.csv
```

---

## Dietary-Stress Cluster Detection

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

---

## Standard Outputs

```text
Phase_II_Host_Harmonized.csv
Phase_II_Microbe_Harmonized.csv
Phase_II_Harmonized_LOR.csv
Phase_II_Cleaned_Signatures.csv
```

---
