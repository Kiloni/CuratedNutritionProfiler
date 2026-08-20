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
