#Install and Load the Package
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("bugsigdbr")

library(bugsigdbr)

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
print(nutrition_conditions)

write.csv(data.frame(Nutrition_Conditions = nutrition_conditions), "Nutrition_Conditions_BugSigDB_Conditions.csv", row.names = FALSE)
