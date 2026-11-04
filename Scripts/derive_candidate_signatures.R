# =============================================================================
# CNP — derive CANDIDATE signatures that fill gaps in Signatures/NutSignaturesList.R.
# The curated list is not modified; candidates are written to
#   Signatures/CNP_Candidate_Signatures.csv   (signature, gene, species, source, method)
# and are validated alongside the curated signatures in
#   Reports/CNP_Signature_Validation.Rmd
#
# Candidates
#   Vitamin_E_Curated_Fixed   curated Vitamin_E list with outdated symbols replaced
#                             (CDC2 -> CDK1, PKC -> PRKCA, the alpha isoform that
#                             alpha-tocopherol inhibits); miRNAs kept but are not
#                             measured by mRNA platforms
#   Vitamin_E_GWAS            protein-coding genes mapped to genome-wide significant
#                             (p <= 5e-8) GWAS Catalog associations for
#                             alpha-tocopherol (EFO_0007898) and tocopherol (EFO_0004738)
#   Vitamin_K_GWAS            same, vitamin K measurement (EFO_0004618); this trait has no
#                             genome-wide significant hits, so the GWAS Catalog's own
#                             inclusion threshold (p <= 1e-5) is used (recorded in `method`)
#   Carnitine_GWAS            same, carnitine measurement (EFO_0010469)
#   Vitamin_K_GSE148011_Up/Down
#                             re-derived from the vitamin K arm of GSE148011 ONLY
#                             (vitamin K precursor 10/20 uM vs 0 uM; one QC-failed library
#                             excluded), limma-trend on log2 CPM: top 25 genes per direction
#                             with FDR < 0.1 (the response is weak: no gene reaches 2-fold).
#                             Intended to replace the curated
#                             Vitamin_K list, which appears to be an annotation artifact.
#
# Also checks that B6_Pyridoxine equals the genes mapped to vitamin B6 measurement
# (EFO_0004621), documenting its source.
#
# Inputs: GWAS trait TSVs from Scripts/download_cnp_source_data.R and
#         data/processed/Vitamin_K/GSE148011_* from Scripts/process_cnp_source_data.R
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(AnnotationDbi); library(org.Hs.eg.db)
})
PROJECT_DIR <- if (file.exists("curatedNutritionDataset.Rproj")) getwd() else
  path.expand("~/RStudio Projects/curatedNutritionDataset")
RAW <- file.path(PROJECT_DIR, "data-raw", "CNP_Source_Data")
GWAS_P <- 5e-8

ex <- parse(file.path(PROJECT_DIR, "Signatures", "NutSignaturesList.R")); eval(ex[[1]])   # NutSignatures

gwas_genes <- function(group, efo_ids) {
  rows <- rbindlist(lapply(efo_ids, function(id) {
    f <- file.path(RAW, group, "GWAS", paste0(id, "_trait_associations.tsv"))
    if (!file.exists(f)) stop("missing ", f, " - run Scripts/download_cnp_source_data.R first")
    x <- fread(f, sep = "\t", quote = "", fill = TRUE)
    x[, .(efo = id, p_value = as.numeric(p_value), mapped_genes)]
  }))
  all_g <- unique(trimws(unlist(strsplit(rows$mapped_genes, "[;,]")))); all_g <- all_g[nzchar(all_g) & !is.na(all_g)]
  pc <- function(r) {
    g <- unique(trimws(unlist(strsplit(r$mapped_genes, "[;,]")))); g <- g[nzchar(g) & !is.na(g)]
    if (!length(g)) return(character(0))
    type <- suppressMessages(mapIds(org.Hs.eg.db, g, "GENETYPE", "SYMBOL", multiVals = "first"))
    sort(g[!is.na(type) & type == "protein-coding"])
  }
  thr <- GWAS_P; sel <- rows[!is.na(p_value) & p_value <= thr]; genes <- pc(sel)
  if (length(genes) < 3) {      # too few genome-wide hits: use the Catalog's own inclusion threshold
    thr <- 1e-5; sel <- rows[!is.na(p_value) & p_value <= thr]; genes <- pc(sel)
  }
  list(genes = genes, n_assoc = nrow(sel), threshold = thr, all = all_g)
}

out <- list()
add <- function(name, genes, species, source, method)
  out[[name]] <<- data.table(signature = name, gene = genes, species = species, source = source, method = method)

# ---- Vitamin E: curated list with symbols corrected -----------------------------
ve <- NutSignatures$Vitamin_E
ve <- unique(sub("^CDC2$", "CDK1", sub("^PKC$", "PRKCA", ve)))
add("Vitamin_E_Curated_Fixed", ve, "human", "NutSignatures$Vitamin_E (curated)",
    "CDC2 -> CDK1, PKC -> PRKCA; MIR122/MIR125B kept (not on mRNA platforms)")

# ---- GWAS-based lists --------------------------------------------------------------
for (spec in list(list("Vitamin_E_GWAS", "Vitamin_E", c("EFO_0007898", "EFO_0004738")),
                  list("Vitamin_K_GWAS", "Vitamin_K", "EFO_0004618"),
                  list("Carnitine_GWAS", "Carnitine", "EFO_0010469"))) {
  g <- gwas_genes(spec[[2]], spec[[3]])
  add(spec[[1]], g$genes, "human", paste("GWAS Catalog", paste(spec[[3]], collapse = " + ")),
      sprintf("protein-coding mapped genes of %d associations with p <= %g", g$n_assoc, g$threshold))
  message(sprintf("%-16s %3d genes from %d associations (p <= %g)", spec[[1]], length(g$genes), g$n_assoc, g$threshold))
}

# ---- B6: confirm provenance of the curated list ---------------------------------
b6 <- gwas_genes("Vitamin_B6", "EFO_0004621")
message("B6_Pyridoxine vs EFO_0004621 mapped genes (any p): identical = ",
        setequal(NutSignatures$B6_Pyridoxine, b6$all), "; in list but not GWAS: ",
        paste(setdiff(NutSignatures$B6_Pyridoxine, b6$all), collapse = ", "),
        "; GWAS but not list: ", paste(setdiff(b6$all, NutSignatures$B6_Pyridoxine), collapse = ", "))

# ---- Vitamin K: re-derive from the vitamin K arm of GSE148011 --------------------
e <- fread(file.path(PROJECT_DIR, "data", "processed", "Vitamin_K", "GSE148011_expr.csv.gz"), data.table = FALSE)
m <- as.matrix(e[, -1]); rownames(m) <- e$gene
p <- read.csv(file.path(PROJECT_DIR, "data", "processed", "Vitamin_K", "GSE148011_pheno.csv"))
p <- p[p$condition %in% c("K3_0", "K3_10", "K3_20"), ]
m <- m[, p$sample]; m <- m[rowMeans(m) > 1, ]                       # expressed genes only
vitk <- factor(ifelse(p$condition == "K3_0", "none", "vitK"), levels = c("none", "vitK"))
fit <- eBayes(lmFit(m, model.matrix(~ vitk)), trend = TRUE)
tt <- topTable(fit, coef = "vitkvitK", number = Inf)
# The response is weak: no gene reaches |log2FC| >= 1 at FDR < 0.05 (max ~1.6-fold), and ~48% of
# genes have p < 0.1, suggesting a broad (possibly technical) shift between arms. Candidates are
# therefore taken at FDR < 0.1 without a fold-change filter and must prove themselves externally.
sig <- tt[tt$adj.P.Val < 0.1, ]
up <- head(rownames(sig)[sig$logFC > 0], 25); dn <- head(rownames(sig)[sig$logFC < 0], 25)
message(sprintf("GSE148011 vitamin K arm: %d genes FDR<0.1 (%d up, %d down); kept %d up, %d down",
                nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0), length(up), length(dn)))
meth <- "limma-trend, vitamin K precursor 10/20 uM (n=7) vs 0 uM (n=4), FDR<0.1 (no gene reaches 2-fold), top 25 by p"
add("Vitamin_K_GSE148011_Up", up, "mouse", "GSE148011 vitamin K arm", paste(meth, "- higher with vitamin K"))
add("Vitamin_K_GSE148011_Down", dn, "mouse", "GSE148011 vitamin K arm", paste(meth, "- lower with vitamin K"))
fwrite(cbind(gene = rownames(tt), tt), file.path(PROJECT_DIR, "Reports", "results", "GSE148011_vitaminK_limma.csv"))

res <- rbindlist(out)
fwrite(res, file.path(PROJECT_DIR, "Signatures", "CNP_Candidate_Signatures.csv"))
print(res[, .(n_genes = .N, example = paste(head(gene, 8), collapse = ", ")), by = signature])
