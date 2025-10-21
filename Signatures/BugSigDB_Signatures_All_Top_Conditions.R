# Install packages for R 4.4.3
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("bugsigdbr", "ComplexHeatmap", "tidyverse", "clusterProfiler"))
install.packages("circlize")

library(bugsigdbr)
library(ComplexHeatmap)
library(tidyverse)
library(clusterProfiler)
library(circlize)
library(grid)

# Import BugSigDB dataset
bsdb <- importBugSigDB()

# Expanded categories list
categories <- list(
  "Metabolic and Endocrine Disorders" = c("Type 1 Diabetes","Type 2 Diabetes","Obesity",
                                          "Metabolic Syndrome","Gestational Diabetes",
                                          "Hyperglycemia","Insulin Resistance"),
  "Gastrointestinal Disorders" = c("Inflammatory Bowel Disease","Crohn’s Disease","Ulcerative Colitis",
                                   "Irritable Bowel Syndrome","Colorectal Cancer","Celiac Disease",
                                   "Gastroesophageal Reflux Disease"),
  "Neurological and Psychiatric Disorders" = c("Parkinson’s Disease","Alzheimer’s Disease",
                                               "Autism Spectrum Disorder","Schizophrenia",
                                               "Depression","Anxiety Disorders","Post-Traumatic Stress Disorder"),
  "Cancer and Tumor Conditions" = c("Breast Cancer","Melanoma","Oral Cancer","Colon Cancer","Prostate Cancer"),
  "Cardiovascular and Circulatory Diseases" = c("Atherosclerosis","Hypertension","Coronary Artery Disease"),
  "Infectious Diseases" = c("COVID-19","HIV/AIDS","Tuberculosis","Helicobacter pylori infection",
                            "Clostridium difficile infection"),
  "Autoimmune and Inflammatory Conditions" = c("Rheumatoid Arthritis","Systemic Lupus Erythematosus",
                                               "Psoriasis","Multiple Sclerosis"),
  "Respiratory Conditions" = c("Asthma","Chronic Obstructive Pulmonary Disease"),
  "Skin and Allergic Conditions" = c("Atopic Dermatitis","Food Allergies","Acne","Eczema"),
  "Reproductive and Hormonal Conditions" = c("Polycystic Ovary Syndrome","Menstrual Cycle Disorders",
                                             "Pregnancy-related complications"),
  "Other Conditions" = c("Kidney Stones","Liver Disease","Sepsis","Cancer Therapy Effects")
)

# Collect enrichment results for summary and combined heatmap
all_enrichment_results <- list()
top_conditions_global <- c()

for (cat_name in names(categories)) {
  conds <- categories[[cat_name]]
  subset_data <- bsdb %>% filter(Condition %in% conds)
  
  if (nrow(subset_data) == 0) next
  
  enrichment_matrix <- subset_data %>%
    select(Condition, `MetaPhlAn taxon names`) %>%
    unnest(cols = `MetaPhlAn taxon names`) %>%
    group_by(`MetaPhlAn taxon names`, Condition) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = count, values_fill = 0)
  
  # Combine low-signal conditions into "Other"
  low_signal <- colSums(enrichment_matrix[,-1]) < 3
  if (any(low_signal)) {
    enrichment_matrix$Other <- rowSums(enrichment_matrix[, names(low_signal)[low_signal]])
    enrichment_matrix <- enrichment_matrix[, c(1, names(low_signal)[!low_signal], "Other")]
  }
  
  enrichment_matrix$total <- rowSums(enrichment_matrix[,-1])
  enrichment_matrix <- enrichment_matrix %>% arrange(desc(total)) %>% slice(1:50)
  
  mat <- as.matrix(enrichment_matrix[,-c(1, ncol(enrichment_matrix))])
  rownames(mat) <- enrichment_matrix$`MetaPhlAn taxon names`
  
  # Log scale for better color contrast
  mat_scaled <- log10(mat + 1)
  
  term2gene <- subset_data %>%
    select(Condition, `MetaPhlAn taxon names`) %>%
    unnest(cols = `MetaPhlAn taxon names`) %>%
    distinct() %>%
    rename(term = Condition, gene = `MetaPhlAn taxon names`)
  
  enrich_res <- enricher(rownames(mat), TERM2GENE = term2gene)
  all_enrichment_results[[cat_name]] <- enrich_res@result
  
  top_conditions_global <- c(top_conditions_global,
                             enrich_res@result %>% arrange(p.adjust) %>% slice(1:3) %>% pull(ID))
  
  sig_terms <- enrich_res@result %>% filter(p.adjust < 0.05) %>% pull(ID)
  
  # Dynamic legend title
  legend_title <- if (length(sig_terms) > 0) {
    "Log10(Enrichment Score)\n(Red = significant taxa)"
  } else {
    "Log10(Enrichment Score)"
  }
  
  # Improved color gradient
  col_fun <- colorRamp2(c(min(mat_scaled), median(mat_scaled), max(mat_scaled)),
                        c("lightblue", "lightyellow", "red"))
  
  row_colors <- ifelse(rownames(mat_scaled) %in% sig_terms, "red", "black")
  
  # Heatmap object with hierarchical clustering and full row names
  ht <- Heatmap(mat_scaled, name = "Log10 Enrichment", col = col_fun,
                cluster_rows = TRUE, cluster_columns = TRUE,
                column_title = cat_name,
                row_names_gp = gpar(fontsize = 8, col = row_colors),
                row_names_max_width = unit(12, "cm"),
                heatmap_legend_param = list(title = legend_title,
                                            legend_direction = "horizontal",
                                            legend_width = unit(4, "cm")))
  
  # Show in RStudio Plot Viewer
  draw(ht, heatmap_legend_side = "bottom")
}

# Combined top conditions heatmap
top_conditions_global <- unique(top_conditions_global)
combined_subset <- bsdb %>% filter(Condition %in% top_conditions_global)

combined_matrix <- combined_subset %>%
  select(Condition, `MetaPhlAn taxon names`) %>%
  unnest(cols = `MetaPhlAn taxon names`) %>%
  group_by(`MetaPhlAn taxon names`, Condition) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = count, values_fill = 0)

combined_matrix$total <- rowSums(combined_matrix[,-1])
combined_matrix <- combined_matrix %>% arrange(desc(total)) %>% slice(1:50)

mat_combined <- as.matrix(combined_matrix[,-c(1, ncol(combined_matrix))])
rownames(mat_combined) <- combined_matrix$`MetaPhlAn taxon names`

mat_combined_scaled <- log10(mat_combined + 1)
col_fun <- colorRamp2(c(min(mat_combined_scaled), median(mat_combined_scaled), max(mat_combined_scaled)),
                      c("lightblue", "lightyellow", "red"))

ht_combined <- Heatmap(mat_combined_scaled, name = "Log10 Enrichment", col = col_fun,
                       cluster_rows = TRUE, cluster_columns = TRUE,
                       column_title = "Top Conditions Across All Categories",
                       row_names_gp = gpar(fontsize = 8),
                       row_names_max_width = unit(20, "cm"),
                       height = unit(nrow(mat_scaled) * 0.18, "cm"), 
                       heatmap_legend_param = list(title = "Log10(Enrichment Score)",
                                                   legend_direction = "horizontal",
                                                   legend_width = unit(3, "cm")))

# Show combined heatmap in RStudio Plot Viewer
draw(ht_combined, heatmap_legend_side = "bottom")

# Export combined enrichment results
combined_enrichment <- bind_rows(all_enrichment_results, .id = "Category")
write.csv(combined_enrichment, "combined_enrichment_results.csv", row.names = FALSE)