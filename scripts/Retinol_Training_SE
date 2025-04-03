#Creating a Summarized experiment with count data and metadata
``
library(TBSignatureProfiler)
library(SummarizedExperiment)
library(tidyverse)
Retinol_SRA_GSE266352_dataset <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/GSE266352_raw_counts_GRCh38.p13_NCBI_GeneSymbol_TRNAVCAC_Combined.csv", header = T) |>
  tibble::column_to_rownames(var = "Symbol")

Retinol_SRA_GSE266352_metadata <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/GSE266352_metadata.csv")
Retinol_metadata <- as.data.frame(Retinol_SRA_GSE266352_metadata[,2:5], row.names = Retinol_SRA_GSE266352_metadata$Sample_ID)
Retinol_SE <- SummarizedExperiment(assays = list(counts = Retinol_SRA_GSE266352_dataset), colData=Retinol_metadata)
``

#Analyzing the file to ensure all metadata is (1) and assays is (counts, log, and counts_to_CPM)
``
dim(Retinol_SE)
assays(Retinol_SE)
Retinol_SE <- mkAssay(Retinol_SE, log = TRUE, counts_to_CPM = TRUE)
``

#Check to make sure the SE looks right
``
Retinol_SE
``

# Save an offline copy to the Nutrition folder on the Mac, just in case!
``
saveRDS(Retinol_SE, file = "Documents/Johnson Lab/Data/Nutrition Data/Retinol_Training_SE.RDS")
``
