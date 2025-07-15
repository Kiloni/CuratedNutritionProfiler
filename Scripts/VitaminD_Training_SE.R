#Creating a Summarized experiment with count data and metadata
``
library(TBSignatureProfiler)
library(SummarizedExperiment)
library(tidyverse)
VitaminD_GSE63086_Data <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/VitaminD_GSE63086_raw_counts_GRCh38.p13_NCBI_Symbol.csv", header = TRUE) |>
  tibble::column_to_rownames(var = "Symbol")

VitaminD_GSE63086_metadata <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/VitaminD_GSE63086_metadata.csv")
VitD_GSE63086_metadata <- as.data.frame(VitaminD_GSE63086_metadata[,2:6], row.names = VitaminD_GSE63086_metadata$Sample_ID)
VitaminD_GSE63086_SE <- SummarizedExperiment(assays = list(counts = VitaminD_Data), colData=VitD_metadata)
``

#Analyzing the file to ensure all metadata is (1) and assays is (counts, log, and counts_to_CPM)
``
dim(VitaminD_GSE63086_SE)
assays(VitaminD_GSE63086_SE)
VitaminD_GSE63086_SE <- mkAssay(VitaminD_GSE63086_SE, log = TRUE, counts_to_CPM = TRUE)
``

#Check to make sure the SE looks right
``
VitaminD_GSE63086_SE
``

# Save an offline copy to the Nutrition folder on the Mac, just in case!
``
saveRDS(VitaminD_GSE63086_SE, file = "Documents/Johnson Lab/Data/Nutrition Data/VitaminD_GSE63086_Training_SE.RDS")
``

#Vitamin D mouse training data: 
#GSE283020_mouse_vitd_def_th2_matrix.csv
#GSE283020_mouse_vitd_metadata.csv


#Creating a Summarized experiment with count data and metadata
``
library(TBSignatureProfiler)
library(SummarizedExperiment)
library(tidyverse)
VitaminD_GSE283020_Mouse_Data <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/GSE283020_mouse_vitd_def_th2_matrix.csv", header = TRUE) |>
  tibble::column_to_rownames(var = "Symbol")

VitaminD_GSE283020_mouse_metadata <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/GSE283020_mouse_vitd_metadata.csv")
VitD_GSE283020_metadata <- as.data.frame(VitaminD_GSE283020_mouse_metadata[,2:4], row.names = VitaminD_GSE283020_mouse_metadata$Sample_ID)
VitaminD_GSE283020_SE <- SummarizedExperiment(assays = list(counts = VitaminD_GSE283020_Mouse_Data), colData=VitD_GSE283020_metadata)
``

#Analyzing the file to ensure all metadata is (1) and assays is (counts, log, and counts_to_CPM)
``
dim(VitaminD_GSE283020_SE)
assays(VitaminD_GSE283020_SE)
VitaminD_GSE283020_SE <- mkAssay(VitaminD_GSE283020_SE, log = TRUE, counts_to_CPM = TRUE)
``

#Check to make sure the SE looks right
``
VitaminD_GSE283020_SE
``

##Creating a boxplot of all vitamin D signatures
# Load required libraries
library(readr)

# Load the data

data <- read_csv("Documents/Johnson Lab/Data/Nutrition Data/VitD mouse training results.csv")

# Rename columns for clarity
colnames(data)[1:3] <- c("Batch", "Genotype", "Treatment")

# Select only Vitamin D signature columns and Treatment
vitd_data <- data %>%
  select(Treatment, starts_with("Vitamin_D_"))

# Reshape the data to long format
vitd_long <- vitd_data %>%
  pivot_longer(cols = starts_with("Vitamin_D_"),
               names_to = "Signature",
               values_to = "Expression")


# Create the boxplot
ggplot(vitd_long, aes(x = Signature, y = Expression, fill = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Differences by Treatment for Vitamin D Signatures",
       x = "Vitamin D Signature",
       y = "Expression")


##Run an ANOVA analysis
# Load necessary libraries
library(tidyverse)
library(car)

# Load the data
data <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/vitamin_d_signatures.csv")

# Convert Treatment to factor
data$Treatment <- as.factor(data$Treatment)

# Perform ANOVA for each Vitamin D signature
vitd_cols <- grep("Vitamin_D", names(data), value = TRUE)

anova_results <- lapply(vitd_cols, function(col) {
  formula <- as.formula(paste(col, "~ Treatment"))
  aov_result <- aov(formula, data = data)
  summary(aov_result)
})
names(anova_results) <- vitd_cols
# Print ANOVA results
anova_results

##Generating AUC/ROC analyses
# Load required libraries
library(pROC)
library(ggplot2)
library(readr)
library(dplyr)

# Load the cleaned dataset
data <- read_csv("Documents/Johnson Lab/Data/Nutrition Data/cleaned_vitamin_d_signatures.csv")

# Create binary classification: 1 = VitD Deficient, 0 = Regular
data <- data %>%
  mutate(VitD_Status = ifelse(Treatment == "VitD Deficient", 1, 0))

# Identify Vitamin D signature columns
vitd_cols <- grep("Vitamin_D", names(data), value = TRUE)

# Loop through each signature and generate ROC curves
for (col in vitd_cols) {
  roc_obj <- roc(data$VitD_Status, data[[col]])
  auc_val <- auc(roc_obj)
  
  # Plot ROC curve
  plot(roc_obj, col = "darkorange", lwd = 2,
       main = paste("ROC Curve -", col),
       legacy.axes = TRUE)
  abline(a = 0, b = 1, col = "navy", lty = 2)
  legend("bottomright", legend = paste("AUC =", round(auc_val, 2)), col = "darkorange", lwd = 2)
}

