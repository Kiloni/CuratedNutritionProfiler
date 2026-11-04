# CuratedNutritionProfiler

**CuratedNutritionProfiler (CNP)** is an R framework for measuring nutrition-related biology in
the transcriptomes of patients with disease. It assembles curated gene-expression signatures
for individual micronutrients, validates each signature against experiments in which nutrient
status was controlled, and then applies the validated signatures to patient cohorts. The aim is
to identify nutrient-linked molecular states that accompany disease and could guide
precision-nutrition and precision-medicine strategies, such as stratifying patients for
targeted supplementation or monitoring their response to it.

The initial disease focus is autoimmune disease (multiple sclerosis, systemic lupus
erythematosus, rheumatoid arthritis) and chronic fatigue syndrome. The framework is also
designed to extend to chronic infection and cancer.

> **Status:** under active development. The signature validation stage is complete for the
> nutrients with suitable public data (see [Validation status](#validation-status)). A
> Bioconductor package and an interactive Shiny interface are planned.

## Approach

CNP works in three stages. A signature moves to the next stage only after it passes the
previous one.

1. **Curate.** Nutrient signatures are assembled from deficiency and supplementation
   experiments, from genome-wide association studies of circulating nutrient levels (GWAS
   Catalog), and from the literature. The curated set is in
   [`Signatures/NutSignaturesList.R`](Signatures/NutSignaturesList.R). Data-derived candidates
   that fill gaps in it are in
   [`Signatures/CNP_Candidate_Signatures.csv`](Signatures/CNP_Candidate_Signatures.csv).
2. **Validate.** Each signature is tested on public studies (GEO) where nutrient status was
   experimentally defined: deficient vs replete diets, or nutrient vs vehicle treatment.
   - **Scoring:** samples are scored with GSVA, ssGSEA, PLAGE, zScore and singscore.
   - **Statistics:** ANOVA, ROC/AUC, logistic regression and random forest.
   - **Internal validation:** leave-one-out cross-validation, plus leave-pair-out AUC, on a
     training dataset.
   - **External validation:** the frozen models are applied to independent datasets from other
     tissues, platforms and species.
   - **Covariates:** known study covariates (genotype, batch, sex) are removed before modelling.
3. **Apply.** Validated signatures are scored in patient transcriptomes and integrated with
   microbiome and dietary data. This identifies disease-associated nutrient states and the
   patients most likely to benefit from a nutritional intervention.

## Validation status

The full results, with figures, are in the knitted report
[`Reports/CNP_Signature_Validation.html`](Reports/CNP_Signature_Validation.html); the underlying
tables are in [`Reports/results/`](Reports/results/). The categories below summarise that
report. "Internal" means the training dataset; "external" means independent datasets.

All 17 curated signatures have now been tested against 47 processed datasets (as of
2026-11-04).

| Status | Signatures |
|---|---|
| Reproduces in independent datasets | `B9_Folate`, `A_Retinol_GLM_43`, `A_Retinol_Quiles_20`, `Carnitine_Mito`, `Iron_Homeostasis`¹ |
| Strong internally, limited external evidence | `B5_Pantothenic Acid`², `Vitamin_C`³, `B6_Pyridoxine`⁴ |
| Weak, inconsistent or algorithm-dependent | `A_Retinol_GWAS`, `A_Retinol_MSigDB_Filtered`, `Vitamin_D_Aigner_35`, `B12_Cobalamin`, `B1_Thiamine`, `B3_Niacin`, `Vitamin_E`, `Vitamin_K`, `Omega_3_FattyAcids` |

¹ Reproduces in animal tissues and cultured cells, but does not separate heart-failure patients
by iron status in whole blood (GSE175739). This is a first indication that tissue signatures
may not transfer directly to blood.
² Confirmed in a single external dataset with three samples per group.
³ Consistent only with the PLAGE score.
⁴ Fails its only external dataset, which uses a B6-pathway inhibitor rather than dietary deficiency.

The report documents the caveats for each signature, including circular training sets, small
groups, treatment-vs-deficiency designs and cross-species mapping. Changes to the curated list
(for example, the 2026-11-04 replacement of `Vitamin_K`) are recorded both in the report and in
the signature file.

## Repository layout

| Path | Contents |
|---|---|
| `Signatures/` | Curated signatures (`NutSignaturesList.R`), data-derived candidates (`CNP_Candidate_Signatures.csv`), microbial signatures (BugSigDB) |
| `Scripts/download_cnp_source_data.R` | Downloads GEO series and GWAS Catalog associations and writes a checksummed manifest |
| `Scripts/process_cnp_source_data.R` | Converts each study into a gene-symbol × sample log2 matrix plus a sample table |
| `Scripts/derive_candidate_signatures.R` | Builds the candidate signatures |
| `Reports/CNP_Signature_Validation.Rmd` | Validation analysis and report (child template `_signature_section.Rmd`) |
| `data/processed/` | Processed expression and phenotype tables (`.csv.gz`), with `dataset_index.csv` |
| `data-raw/CNP_Source_Data/` | Download manifest, log and GWAS tables. Raw GEO files are deleted after processing and can be re-fetched. |
| `Skills/` | Phase-by-phase specifications for the broader multi-omics framework |

## Reproducing the analysis

Requirements: R ≥ 4.1 with Bioconductor packages GSVA, singscore, limma, edgeR,
AnnotationDbi, org.Hs.eg.db, org.Mm.eg.db and hugene20sttranscriptcluster.db, plus CRAN
packages data.table, readxl, pROC, randomForest, babelgene, ggplot2, kableExtra and rmarkdown.

```r
source("Scripts/download_cnp_source_data.R")     # fetch raw data (skips processed and excluded accessions)
source("Scripts/process_cnp_source_data.R")      # build data/processed/
source("Scripts/derive_candidate_signatures.R")  # build candidate signatures
rmarkdown::render("Reports/CNP_Signature_Validation.Rmd")
```

Every download is listed with its source URL and MD5 checksum in
`data-raw/CNP_Source_Data/download_manifest.csv`, so processed data can be traced to, and
re-created from, the original public records.

## Data sources

Gene expression data are from the NCBI Gene Expression Omnibus (GEO). Genetic associations are
from the NHGRI-EBI GWAS Catalog. Accession numbers, study designs and the case/control
definitions used for each study are listed in the validation report.

## Contact

Developed in the Johnson Lab, Rutgers New Jersey Medical School. For questions, please open an
issue on this repository.
