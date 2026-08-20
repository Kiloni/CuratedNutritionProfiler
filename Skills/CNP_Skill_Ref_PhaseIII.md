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
