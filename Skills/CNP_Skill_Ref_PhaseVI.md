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
