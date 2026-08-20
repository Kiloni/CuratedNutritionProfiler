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
