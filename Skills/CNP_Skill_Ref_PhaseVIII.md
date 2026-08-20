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
---
