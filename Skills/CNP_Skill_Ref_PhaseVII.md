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
