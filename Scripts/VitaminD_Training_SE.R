#Creating a Summarized experiment with count data and metadata
``
library(TBSignatureProfiler)
library(SummarizedExperiment)
library(tidyverse)
VitaminD_Data <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/VitaminD_GSE63086_raw_counts_GRCh38.p13_NCBI_Symbol.csv", header = TRUE) |>
  tibble::column_to_rownames(var = "Symbol")

VitaminD_GSE63086_metadata <- read.csv("Documents/Johnson Lab/Data/Nutrition Data/VitaminD_GSE63086_metadata.csv")
VitD_metadata <- as.data.frame(VitaminD_GSE63086_metadata[,2:6], row.names = VitaminD_GSE63086_metadata$Sample_ID)
VitaminD_SE <- SummarizedExperiment(assays = list(counts = VitaminD_Data), colData=VitD_metadata)
``

#Analyzing the file to ensure all metadata is (1) and assays is (counts, log, and counts_to_CPM)
``
dim(VitaminD_SE)
assays(VitaminD_SE)
VitaminD_SE <- mkAssay(VitaminD_SE, log = TRUE, counts_to_CPM = TRUE)
``

#Check to make sure the SE looks right
``
VitaminD_SE
``

# Save an offline copy to the Nutrition folder on the Mac, just in case!
``
saveRDS(VitaminD_SE, file = "Documents/Johnson Lab/Data/Nutrition Data/VitaminD_Training_SE.RDS")
``
