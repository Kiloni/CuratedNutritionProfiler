#Install and Load the Package
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("bugsigdbr")

library(bugsigdbr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(scales)

#Import the BugSigDB Dataset
bsdb <- importBugSigDB()

#Extract Unique Conditions
all_conditions <- unique(bsdb$Condition)


# Define nutrition-related keywords
nutrition_keywords <- c("diet", "nutrition", "fiber", "fat", "carbohydrate", "protein", "fasting", "calorie", "polyphenol", "vitamin", "supplement")

# Filter conditions using keyword matching (case-insensitive)
nutrition_conditions <- all_conditions[
  grepl(paste(nutrition_keywords, collapse = "|"), all_conditions, ignore.case = TRUE)
]

# Print and save a CSV of the nutrition-related conditions                                 
print(nutrition_conditions)
write.csv(data.frame(Nutrition_Conditions = nutrition_conditions), "data-raw/Microbes/Nutrition_Conditions_BugSigDB_Conditions.csv", row.names = FALSE)

# Filter studies with nutrition-related conditions
nutrition_studies <- bsdb[grepl(paste(nutrition_keywords, collapse = "|"), bsdb$Condition, ignore.case = TRUE),
]

# Extract microbial signatures (species-level)
nutrition_signatures <- getSignatures(nutrition_studies, tax.id.type = "taxname", tax.level = "species")

# View the signature results
str(nutrition_signatures)

# Save to a GMT file for enrichment analysis
writeGMT(nutrition_signatures, "Signatures/nutrition_signatures_bugsigdb.gmt")

# Define microbial signatures for each condition
nutrition_signatures <- list(
  Non_alcoholic_fatty_liver_disease = c("Bacteroides fragilis", "Faecalibacterium prausnitzii"),
  Chronic_fatigue_syndrome = c("Bacteroides fragilis", "Alistipes putredinis", "Faecalibacterium prausnitzii"),
  High_fat_diet = c("Bilophila wadsworthia", "Alistipes putredinis", "Bacteroides fragilis"),
  Diet = c("Roseburia intestinalis", "Bifidobacterium adolescentis"),
  Response_to_diet = c("Faecalibacterium prausnitzii", "Eubacterium rectale"),
  Body_fat_percentage = c("Akkermansia muciniphila", "Lactobacillus plantarum"),
  Response_to_ketogenic_diet = c("Bacteroides fragilis", "Roseburia intestinalis"),
  Fasting = c("Faecalibacterium prausnitzii", "Bifidobacterium adolescentis"),
  Diet_measurement = c("Lactobacillus plantarum", "Eubacterium rectale"),
  Response_to_supplemental_oxygen = c("Alistipes putredinis", "Bilophila wadsworthia"),
  Diarrhea = c("Bacteroides fragilis", "Lactobacillus plantarum"),
  Infantile = c("Bifidobacterium adolescentis", "Roseburia intestinalis"),
  Malnutrition = c("Faecalibacterium prausnitzii", "Eubacterium rectale")
)

# Create full microbe-condition presence matrix
all_microbes <- sort(unique(unlist(nutrition_signatures)))
all_conditions <- names(nutrition_signatures)

presence_matrix <- sapply(all_conditions, function(cond) {
  as.integer(all_microbes %in% nutrition_signatures[[cond]])
})
rownames(presence_matrix) <- all_microbes

# View the matrix
print(presence_matrix)

# Optional: Save to CSV
write.csv(presence_matrix, "data-raw/Microbes/BugSigDB_microbe_presence_matrix.csv", row.names = TRUE)

# Perform Fisher's Exact Test
pvals <- matrix(NA, nrow = nrow(presence_matrix), ncol = ncol(presence_matrix),
                dimnames = list(rownames(presence_matrix), colnames(presence_matrix)))

for (i in 1:nrow(presence_matrix)) {
  for (j in 1:ncol(presence_matrix)) {
    a <- presence_matrix[i, j]
    b <- sum(presence_matrix[i, ]) - a
    c <- sum(presence_matrix[, j]) - a
    d <- sum(presence_matrix) - (a + b + c)
    test <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    pvals[i, j] <- test$p.value
  }
}

# Replace 0s and transform to -log10(p-values)
pvals[pvals == 0] <- .Machine$double.xmin
log_pvals <- -log10(pvals)

# Convert matrix to long format
log_pvals_long <- as.data.frame(as.table(log_pvals))
colnames(log_pvals_long) <- c("Microbe", "Condition", "logP")

# Normalize values for color scale
log_pvals_long$fill_color <- rescale(log_pvals_long$logP, to = c(0, 1))

# Determine text color based on fill intensity
log_pvals_long$text_color <- ifelse(log_pvals_long$fill_color > 0.5, "black", "white")

# Plot heatmap
ggplot(log_pvals_long, aes(x = Condition, y = Microbe, fill = logP)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = round(logP, 2), color = text_color), size = 3) +
  scale_fill_gradientn(colors = c("navy", "skyblue", "red"), name = "-log10(p-value)") +
  scale_color_identity() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(title = "-log10(P-values) from Fisher's Exact Test (Nutrition Conditions)",
       x = "Condition", y = "Microbes")

#_________________

# Repeat the process above with all clinical and disease state data
clinical_conditions <- c(
  "Type 1 Diabetes", "Type 2 Diabetes", "Obesity", "Metabolic Syndrome", "Gestational Diabetes",
  "Hyperglycemia", "Insulin Resistance", "Inflammatory Bowel Disease", "Crohn’s Disease",
  "Ulcerative Colitis", "Irritable Bowel Syndrome", "Colorectal Cancer", "Celiac Disease",
  "Gastroesophageal Reflux Disease", "Parkinson’s Disease", "Alzheimer’s Disease",
  "Autism Spectrum Disorder", "Schizophrenia", "Depression", "Anxiety Disorders",
  "Post-Traumatic Stress Disorder", "Breast Cancer", "Melanoma", "Oral Cancer", "Colon Cancer",
  "Prostate Cancer", "Atherosclerosis", "Hypertension", "Coronary Artery Disease", "COVID-19",
  "HIV/AIDS", "Tuberculosis", "Helicobacter pylori infection", "Clostridium difficile infection",
  "Rheumatoid Arthritis", "Systemic Lupus Erythematosus", "Psoriasis", "Multiple Sclerosis",
  "Asthma", "Chronic Obstructive Pulmonary Disease", "Atopic Dermatitis", "Food Allergies",
  "Acne", "Eczema", "Polycystic Ovary Syndrome", "Menstrual Cycle Disorders",
  "Pregnancy-related complications", "Kidney Stones", "Liver Disease", "Sepsis",
  "Cancer Therapy Effects"
)

# Filter studies with clinical/disease conditions
clinical_studies <- bsdb[bsdb$Condition %in% clinical_conditions, ]

# Extract microbial signatures (species-level)
clinical_signatures <- getSignatures(clinical_studies, tax.id.type = "taxname", tax.level = "species")

# Save to a GMT file for enrichment analysis
writeGMT(clinical_signatures, "Signatures/clinical_signatures_bugsigdb.gmt")

# Filter BugSigDB for relevant conditions
filtered_bsdb <- bsdb %>% filter(Condition %in% clinical_conditions)

# Extract species-level signatures
clinical_signatures_filtered <- getSignatures(filtered_bsdb, tax.id.type = "taxname", tax.level = "species")

# Create presence matrix
all_microbes <- sort(unique(unlist(clinical_signatures_filtered)))
presence_matrix <- sapply(names(clinical_signatures_filtered), function(cond) {
  as.integer(all_microbes %in% clinical_signatures_filtered[[cond]])
})
rownames(presence_matrix) <- all_microbes
colnames(presence_matrix) <- names(clinical_signatures_filtered)

# Perform Fisher's Exact Test
pvals <- matrix(NA, nrow = nrow(presence_matrix), ncol = ncol(presence_matrix),
                dimnames = list(rownames(presence_matrix), colnames(presence_matrix)))

for (i in 1:nrow(presence_matrix)) {
  for (j in 1:ncol(presence_matrix)) {
    a <- presence_matrix[i, j]
    b <- sum(presence_matrix[i, ]) - a
    c <- sum(presence_matrix[, j]) - a
    d <- sum(presence_matrix) - (a + b + c)
    test <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    pvals[i, j] <- test$p.value
  }
}

# Convert to -log10(p-values)

pvals[pvals == 0] <- .Machine$double.xmin
log_pvals <- -log10(pvals)

# Plot heatmap
pheatmap(log_pvals,
         color = colorRampPalette(c("navy", "skyblue", "red"))(100),
         main = "-log10(P-values) from Fisher's Exact Test (Clinical Conditions)",
         fontsize_row = 6,
         fontsize_col = 8,
         angle_col = 45)

# Plot heatmap
pheatmap(log_pvals,
         color = colorRampPalette(c("navy", "skyblue", "red"))(100),
         main = "-log10(P-values) from Fisher's Exact Test (Clinical Conditions)",
         fontsize_row = 6,               # Smaller font for many microbes
         fontsize_col = 8,               # Smaller font for many conditions
         angle_col = 45,                 # Rotate column labels
         display_numbers = FALSE,        # Hide numbers to reduce clutter
         border_color = "grey80",
         legend = TRUE,
         treeheight_row = 0,             # Disable row clustering tree
         treeheight_col = 0)             # Disable column clustering tree
