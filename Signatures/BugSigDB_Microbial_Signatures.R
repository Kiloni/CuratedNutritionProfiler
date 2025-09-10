#Install and Load the Package
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("bugsigdbr")

library(bugsigdbr)
library(tidyverse)
library(reshape2)
library(pheatmap)

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
write.csv(data.frame(Nutrition_Conditions = nutrition_conditions), "Nutrition_Conditions_BugSigDB_Conditions.csv", row.names = FALSE)

# Filter studies with nutrition-related conditions
nutrition_studies <- bsdb[grepl(paste(nutrition_keywords, collapse = "|"), bsdb$Condition, ignore.case = TRUE),
]

# Extract microbial signatures (species-level)
nutrition_signatures <- getSignatures(nutrition_studies, tax.id.type = "taxname", tax.level = "species")

# View the signature results
str(nutrition_signatures)

# Save to a GMT file for enrichment analysis
writeGMT(nutrition_signatures, "nutrition_signatures_bugsigdb.gmt")

