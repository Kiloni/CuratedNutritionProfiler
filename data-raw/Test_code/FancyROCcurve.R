# Load required libraries
library(readxl)
library(tidyverse)
library(pROC)
library(dplyr)
library(writexl)

# Load the Excel file with the nutritional dataset
VitCdf <- read_excel("~/Downloads/GSE233598_tpm.xlsx", sheet = "tpm")

# Combine duplicate rows based on the 'ENSEMBL' column by summing their values
VitCdf_combined <- VitCdf %>%
  group_by(ENSEMBL) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

# Remove rows where all read values (excluding 'ENSEMBL' and 'SYMBOL') are zero
read_columns <- setdiff(names(VitCdf_combined), c("ENSEMBL", "SYMBOL"))
VitCdf_cleaned <- VitCdf_combined %>%
  filter(rowSums(select(., all_of(read_columns)) != 0) > 0)

# Save the cleaned dataset to a new Excel file
write_xlsx(df_cleaned, "GSE233598_tpm_cleaned.xlsx")


# Load TPM expression data
tpm_data <- read_excel("~/Downloads/GSE233598_tpm_cleaned.xlsx", sheet = "tpm")

# Clean and format expression data
tpm_data <- tpm_data %>% filter(!is.na(SYMBOL)) %>%
  column_to_rownames("SYMBOL") %>%
  select(-ENSEMBL)

expression_matrix <- as.data.frame(t(tpm_data))

# Load metadata from series matrix
metadata_lines <- readLines("~/Downloads/GSE233598_series_matrix.txt")

# Extract sample titles and treatments
sample_titles <- metadata_lines[grepl("!Sample_title", metadata_lines)] %>%
  str_split("\t") %>% unlist() %>% .[-1]

treatment_lines <- metadata_lines[grepl("!Sample_characteristics_ch1", metadata_lines) & grepl("treatment", metadata_lines)]
treatments <- treatment_lines %>%
  str_split("\t") %>% unlist() %>% .[-1] %>%
  str_replace("treatment: ", "")

# Create metadata dataframe
metadata <- data.frame(
  Sample = rownames(expression_matrix),
  Title = sample_titles,
  Treatment = treatments,
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$Sample

# Define vitamin C signatures
vitamin_c_signatures <- c("SLC23A1", "SLC23A3", "CHPT1", "BCAS3", "SNRPF", "RER1", "MAF", "GSTA5", 
                          "RGS14", "AKT1", "FADS1")
vitamin_c_signatures <- vitamin_c_signatures[vitamin_c_signatures %in% colnames(expression_matrix)]

# Filter samples for ROC analysis (0.40% vitC = 1, no vitC = 0)
metadata$DietBinary <- case_when(
  str_detect(metadata$Treatment, "0.40%") ~ 1,
  str_detect(metadata$Treatment, "no vitC") ~ 0,
  TRUE ~ NA_real_
)

valid_samples <- metadata %>% filter(!is.na(DietBinary)) %>% rownames()

# Prepare data for ROC
roc_data <- lapply(vitamin_c_signatures, function(gene) {
  predictor <- expression_matrix[valid_samples, gene]
  response <- metadata[valid_samples, "DietBinary"]
  
  roc_obj <- roc(response, predictor)
  data.frame(
    fpr = rev(roc_obj$specificities),
    tpr = rev(roc_obj$sensitivities),
    Signature = paste0(gene, " (AUC = ", round(auc(roc_obj), 2), ")")
  )
}) %>% bind_rows()

# Plot ROC curves
ggplot(roc_data, aes(x = fpr, y = tpr, color = Signature)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curves for Vitamin C Signatures",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Signature"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
