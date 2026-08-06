# =====================================================================
# STEP 1: Install Required Packages
# =====================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install the database package and data manipulation tools
# (Press "a" to update all packages if prompted in the console)
BiocManager::install(c("curatedMetagenomicData", "dplyr", "TreeSummarizedExperiment"), ask = FALSE)

library(curatedMetagenomicData)
library(dplyr)
library(TreeSummarizedExperiment)

# =====================================================================
# STEP 2 & 3: Define EXACT Study Names and Download
# =====================================================================
# Corrected the name for the IBDMDB dataset!
target_studies <- c("HMP_2019_ibdmdb", "FengQ_2015")

# Append '.relative_abundance' to tell the database exactly what data type we want
dataset_queries <- paste0(target_studies, ".relative_abundance")

# Collapse them into a single search pattern
search_pattern <- paste(dataset_queries, collapse = "|")

message("Fetching datasets directly from ExperimentHub...")
message("Query: ", search_pattern)

# Download the data
processed_datasets <- curatedMetagenomicData(
  pattern = search_pattern, 
  dryrun = FALSE, 
  counts = FALSE
)

message("\nDownload complete! ", length(processed_datasets), " dataset(s) retrieved.")

# =====================================================================
# STEP 4: Extract Matrices
# =====================================================================
for (study in names(processed_datasets)) {
  message("\nExtracting matrices for: ", study)
  
  # 1. Get the taxonomic abundance matrix
  abundance_matrix <- assay(processed_datasets[[study]])
  
  # 2. Get the clinical/patient metadata
  patient_metadata <- as.data.frame(colData(processed_datasets[[study]]))
  
  # Save them locally as CSV files
  write.csv(abundance_matrix, paste0(study, "_taxonomic_abundance.csv"))
  write.csv(patient_metadata, paste0(study, "_patient_metadata.csv"))
  
  message("Saved: ", paste0(study, "_taxonomic_abundance.csv"))
}

# =====================================================================
# STEP 4: Extract Data into Standard DataFrames (Optional)
# =====================================================================
# The data downloads as 'TreeSummarizedExperiment' objects (the Bioconductor standard).
# Here is how you extract the combined count matrix and the patient metadata into normal R dataframes.

for (study in names(processed_datasets)) {
  message("\nExtracting matrices for: ", study)
  
  # 1. Get the taxonomic abundance matrix (Rows = Microbes, Columns = Samples)
  abundance_matrix <- assay(processed_datasets[[study]])
  
  # 2. Get the clinical/patient metadata (Rows = Samples, Columns = Variables like age, disease, etc.)
  patient_metadata <- as.data.frame(colData(processed_datasets[[study]]))
  
  # Save them locally as CSV files so you have the single combined files you wanted
  write.csv(abundance_matrix, paste0(study, "_taxonomic_abundance.csv"))
  write.csv(patient_metadata, paste0(study, "_patient_metadata.csv"))
  
  message("Saved: ", paste0(study, "_taxonomic_abundance.csv"))
}

message("\nAll available processed data successfully downloaded and combined!")
