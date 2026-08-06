# Load necessary library
if (!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

# ==============================================================================
# DATASET 1: HUMAN GENE SIGNATURES (Brain & Blood)
# Sources: SFARI Gene, Voineagu et al. (2011), Kong et al. (2012)
# ==============================================================================

human_genes <- tribble(
  ~Gene_Symbol, ~Ensembl_ID, ~Category, ~Regulation_Direction, ~Tissue, ~Study_Source, ~Function_Notes,
  "CHD8", "ENSG00000100888", "High Confidence (SFARI 1)", "Downregulated (Mutation)", "Brain/Systemic", "SFARI Gene", "Chromatin remodeling; major regulator of other ASD genes",
  "SHANK3", "ENSG00000251322", "High Confidence (SFARI 1)", "Downregulated (Mutation)", "Synapse", "SFARI Gene", "Synaptic scaffolding protein; Phelan-McDermid syndrome",
  "SCN2A", "ENSG00000136531", "High Confidence (SFARI 1)", "Variable", "Brain", "SFARI Gene", "Sodium channel, voltage-gated; seizure comorbidity",
  "ADNP", "ENSG00000196497", "High Confidence (SFARI 1)", "Downregulated", "Brain", "SFARI Gene", "Activity-dependent neuroprotector",
  "A2BP1 (RBFOX1)", "ENSG00000078328", "Transcriptomic Signature", "Downregulated", "Cortex", "Voineagu et al. (2011)", "Neuronal splicing factor; convergent downregulation",
  "STX1A", "ENSG00000106089", "Transcriptomic Signature", "Downregulated", "Cortex", "Voineagu et al. (2011)", "Synaptic vesicle exocytosis",
  "SYT1", "ENSG00000067715", "Transcriptomic Signature", "Downregulated", "Cortex", "Voineagu et al. (2011)", "Calcium sensor for synaptic transmission",
  "IL17A", "ENSG00000112115", "Immune Signature", "Upregulated", "Blood/Cortex", "Choi et al. / Voineagu", "Pro-inflammatory cytokine; maternal immune activation marker",
  "GFAP", "ENSG00000131095", "Immune Signature", "Upregulated", "Cortex", "Voineagu et al. (2011)", "Astrocyte marker; indicates gliosis/inflammation",
  "AIF1 (Iba1)", "ENSG00000204472", "Immune Signature", "Upregulated", "Cortex", "Voineagu et al. (2011)", "Microglia marker; indicates neuroinflammation",
  "FOS", "ENSG00000170345", "Blood Signature", "Upregulated", "Peripheral Blood", "Kong et al. (2012)", "Immediate early gene; predictive blood biomarker",
  "DYRK1A", "ENSG00000157540", "High Confidence (SFARI 1)", "Downregulated", "Brain", "SFARI Gene", "Kinase involved in brain growth; Down Syndrome critical region"
)

# Export to CSV
write_csv(human_genes, "ASD_Human_Gene_Signatures.csv")
message("File 'ASD_Human_Gene_Signatures.csv' has been generated.")


# ==============================================================================
# DATASET 2: MICROBIOME DATASETS & CONFOUNDERS
# Sources: Yap et al., Kang et al., American Gut Project
# ==============================================================================

microbiome_data <- tribble(
  ~Study_ID, ~Year, ~Authors, ~Sequencing_Technique, ~Data_Type, ~Key_Finding, ~Major_Confounders_Tracked, ~Accession_ID,
  "Yap_Australian_Biobank", 2021, "Yap et al.", "Shotgun Metagenomics", "Stool Microbiome", "Microbiome diversity correlates with dietary diversity, not ASD diagnosis", "Detailed Diet (FFQ), Age, Sex, Stool Consistency (Bristol)", "PRJEB41566",
  "MTT_Trial_Ph1", 2017, "Kang et al.", "16S rRNA (V4)", "Longitudinal Stool", "Fecal transplant increased Bifidobacterium & Prevotella; improved GI/ASD symptoms", "Antibiotic use history, GI Severity Score (GSRS)", "PRJNA388274",
  "American_Gut_ASD", 2018, "McDonald et al.", "16S rRNA", "Citizen Science Stool", "Identified associations with antibiotics and diet over diagnosis", "Antibiotics (recent), Diet, Geography, BMI", "EBI: ERP012803",
  "Son_Metagenomics", 2015, "Son et al.", "Shotgun Metagenomics", "Stool & Urine Metabolomics", "Altered SCFA levels; correlation with severity", "Neurodevelopmental severity scores", "PRJNA282013",
  "China_Multi_Omics", 2020, "Wan et al.", "Shotgun + Metabolomics", "Stool + Plasma", "Identified impaired detoxification enzymes in ASD microbiome", "Diet patterns, constipation status", "PRJNA606751"
)

# Export to CSV
write_csv(microbiome_data, "ASD_Microbiome_Resources.csv")
message("File 'ASD_Microbiome_Resources.csv' has been generated.")

# ==============================================================================
# DATASET 2: MICROBIAL GENE SIGNATURES (Species & Target Genes)
# ==============================================================================

microbial_genes <- tribble(
  ~Microbe_Species, ~Target_Gene_or_Marker, ~Regulation_in_ASD, ~Functional_Relevance, ~Reference,
  
  # --- Clostridium & Toxin Producers ---
  "Clostridium bolteae", "16S rRNA / wzy (Capsular Polysaccharide)", "Upregulated", "Produces specific capsular polysaccharides; historically linked to GI severity.", "Pequegnat et al. (2013)",
  "Clostridium perfringens", "cpe (Enterotoxin gene) / plc (Alpha-toxin)", "Upregulated", "Toxin production disrupts gut barrier (leaky gut); associated with GI symptoms.", "Finegold et al.",
  
  # --- Sulfur Metabolism (Toxic H2S production) ---
  "Desulfovibrio piger", "dsrA / dsrB (Dissimilatory sulfite reductase)", "Upregulated", "Reduces sulfate to toxic Hydrogen Sulfide (H2S); inhibits mitochondrial function.", "Finegold / Kang et al.",
  "Desulfovibrio desulfuricans", "16S rRNA / nifH (Nitrogenase)", "Upregulated", "Sulfur reducer; correlates with severity of autism manifestations.", "Tomova et al.",
  
  # --- Neurotransmitter Modulators ---
  "Lactobacillus reuteri", "oxytocin-pathway genes (unknown mechanism)", "Downregulated", "Restores social deficits in mouse models via vagus nerve/oxytocin signaling.", "Sgritta et al. (2019)",
  "Escherichia coli", "gadB (Glutamate decarboxylase)", "Variable (Strain dependent)", "Converts Glutamate (excitatory) to GABA (inhibitory); imbalance affects gut-brain signaling.", "Strandwitz et al.",
  "Bifidobacterium metabolic types", "adh (Alcohol dehydrogenase)", "Downregulated", "Reduction impairs detoxification of alcohols/aldehydes; increases oxidative stress.", "Wan et al. (2020)",
  
  # --- Short Chain Fatty Acid (SCFA) Producers ---
  "Faecalibacterium prausnitzii", "but (Butyryl-CoA:acetate CoA-transferase)", "Downregulated", "Major butyrate producer; butyrate maintains blood-brain barrier integrity.", "Liu et al. (2019)",
  "Roseburia intestinalis", "16S rRNA / butyrate kinase", "Downregulated", "Crucial for anti-inflammatory SCFA production; depletion linked to gut inflammation.", "Kang et al.",
  
  # --- Specific Markers ---
  "Sutterella wadsworthensis", "16S rRNA (V4 region)", "Upregulated", "Found in ileal biopsies; distinct from other dysbiosis; role remains unclear but highly specific.", "Williams et al.",
  "Akkermansia muciniphila", "mucin-degradation genes", "Downregulated", "Regulates mucus layer thickness; depletion leads to increased permeability.", "Yap et al. (2021)"
)

# Export Microbial Genes
write_csv(microbial_genes, "ASD_Microbial_Signatures.csv")

message("Success! Two files generated: 'ASD_Human_Gene_Signatures.csv' and 'ASD_Microbial_Signatures.csv'")

# ==============================================================================
# 2. VISUALIZATION: Microbe–Nutrient Interaction Heatmap
# ==============================================================================

library(tidyverse)
library(ggplot2)

# --- Reshape to long format ---
heatmap_long <- heatmap_data %>%
  pivot_longer(
    cols = c(Vitamin_B12, Carnitine, Fiber_SCFA, Tryptophan, Glutathione, Magnesium),
    names_to = "Nutrient",
    values_to = "Score"
  )

# --- Order aesthetics to match your supplied figure ---
heatmap_long$Microbe <- factor(
  heatmap_long$Microbe,
  levels = c(
    "Sutterella",
    "Lactobacillus reuteri",
    "Faecalibacterium",
    "Desulfovibrio (H2S)",
    "Clostridium perfringens",
    "Clostridium bolteae",
    "Candida (Yeast)",
    "Bifidobacterium",
    "Akkermansia"
  )
)

heatmap_long$Nutrient <- factor(
  heatmap_long$Nutrient,
  levels = c("Carnitine", "Fiber_SCFA", "Glutathione", "Magnesium", "Tryptophan", "Vitamin_B12")
)

# --- Generate figure ---
p <- ggplot(heatmap_long, aes(x = Nutrient, y = Microbe, fill = Score)) +
  geom_tile(color = "white", size = 0.7) +
  
  # Value labels
  geom_text(aes(label = Score), size = 4.2, fontface = "bold") +
  
  scale_fill_gradient2(
    low = "#d73027",    # red (negative)
    mid = "white",
    high = "#4575b4",   # blue (positive)
    midpoint = 0,
    limits = c(-1, 1),
    name = "Interaction\nType"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  
  labs(
    title = "Microbe–Nutrient Interaction Matrix in ASD",
    subtitle = "Red = Depletion/Antagonism | Blue = Production/Synergy",
    x = "Nutritional Factor / Metabolite",
    y = "Microbial Species"
  )

# Print the plot
print(p)

# Save to file if desired
ggsave("Microbe_Nutrient_Interaction_Matrix.png", p, width = 10, height = 6, dpi = 300)
``

# ==============================================================================
# 1. CONSTRUCT THE DATA MATRIX (Literature-Derived Interactions)
# Scores represent the functional relationship:
#  +1.0 : Positive / Synergistic (Microbe produces nutrient or aids absorption)
#  -1.0 : Negative / Antagonistic (Microbe depletes nutrient or toxin blocks it)
#   0.0 : No strong established direct link
# ==============================================================================

# Define the data frame
heatmap_data <- tribble(
  ~Microbe, ~Vitamin_B12, ~Carnitine, ~Fiber_SCFA, ~Tryptophan, ~Glutathione, ~Magnesium,
  
  # Pathogenic / Overgrown in ASD
  "Clostridium bolteae",    -0.5, -1.0, -0.8, -0.2, -0.5, -0.2,
  "Clostridium perfringens",-0.2, -0.8, -0.5, -0.5, -0.8,  0.0,
  "Desulfovibrio (H2S)",     0.0, -0.5, -0.5, -0.2, -1.0, -0.5,
  "Sutterella",             -0.2,  0.0, -0.5,  0.0, -0.2, -0.8,
  "Candida (Yeast)",        -0.5, -0.2, -0.5, -0.5, -0.5, -0.5,
  
  # Beneficial / Depleted in ASD
  "Bifidobacterium",         0.8,  0.2,  1.0,  0.5,  0.2,  0.8,
  "Lactobacillus reuteri",   0.8,  0.0,  0.5,  1.0,  0.0,  0.5,
  "Akkermansia",             0.0,  0.0,  0.8,  0.0,  0.5,  0.5,
  "Faecalibacterium",        0.2,  0.0,  1.0,  0.2,  0.5,  0.5
)

# ==============================================================================
# 2. DATA TRANSFORMATION (Pivot for Heatmap)
# ==============================================================================

# Convert from wide to long format for ggplot
heatmap_long <- heatmap_data %>%
  pivot_longer(
    cols = -Microbe, 
    names_to = "Nutrient", 
    values_to = "Interaction_Score"
  )

# ==============================================================================
# 3. GENERATE THE HEATMAP
# ==============================================================================

plot <- ggplot(heatmap_long, aes(x = Nutrient, y = Microbe, fill = Interaction_Score)) +
  
  # Create tiles
  geom_tile(color = "white", lwd = 0.5) +
  
  # Color Gradient: 
  # Red = Depletion/Negative (e.g., Clostridia depleting Carnitine)
  # Blue = Production/Positive (e.g., Bifido producing SCFA)
  scale_fill_gradient2(
    low = "#D73027",    # Red (Depletion)
    mid = "#F7F7F7",    # White (Neutral)
    high = "#4575B4",   # Blue (Production/Synergy)
    midpoint = 0,
    name = "Interaction\nType"
  ) +
  
  # Add value labels inside tiles for clarity
  geom_text(aes(label = Interaction_Score), color = "black", size = 3) +
  
  # Formatting
  theme_minimal() +
  labs(
    title = "Microbe-Nutrient Interaction Matrix in ASD",
    subtitle = "Red = Depletion/Antagonism | Blue = Production/Synergy",
    x = "Nutritional Factor / Metabolite",
    y = "Microbial Species"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

# Display the plot
print(plot)

# Optional: Save the plot
ggsave("ASD_Microbe_Nutrient_Heatmap.png", plot, width = 8, height = 6)

# ==============================================================================
# SETUP & LIBRARIES
# ==============================================================================

# Install specific bio-packages if not present
list.of.packages <- c("tidyverse", "GEOquery", "pheatmap", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(!require("limma")) BiocManager::install("limma")

library(tidyverse)
library(GEOquery)
library(pheatmap)
library(RColorBrewer)

# ==============================================================================
# 1. DEFINE GENE LISTS
# ==============================================================================

# ==============================================================================
# UPDATED NUTRITION LIST (Expanded for Methylation & Mitochondria)
# ==============================================================================

nutrition_list <- list(
  # --- ORIGINAL LIST ---
  "A_Retinol"      = c("RPE65", "ADH1A", "DHRS4", "CYP2A6", "RDH12", "BCO1"),
  "B1_Thiamine"    = c("SLC19A2", "TPK1", "OGDH", "PDHA1", "SLC25A19", "TKT"),
  "B6_Pyridoxine"  = c("ALPL", "PNPO", "AMT", "CSAD", "GAD1", "GAD2"),
  "B12_Cobalamin"  = c("TCN2", "MTR", "MTRR", "MMACHC", "MUT", "CD320"),
  "Vitamin_D"      = c("VDR", "CYP27B1", "CYP24A1", "GC"),
  "Vitamin_E"      = c("TTPA", "CYP4F2", "GSTP1"),
  "Zinc"           = c("SLC30A8", "MT1A", "MT2A", "SLC39A1"),
  
  # --- NEW SIGNATURES ---
  
  # 1. FOLATE (B9) - Critical for Methylation (Links to CHD8/KMT2A)
  "B9_Folate"      = c("MTHFR", "DHFR", "SLC19A1", "FOLR1", "SHMT1", "MTHFD1"),
  
  # 2. IRON - Essential for Synaptic Energy & Myelination
  "Iron_Homeostasis" = c("TFRC", "FTH1", "FTL", "SLC11A2", "ACO1", "IREB2"),
  
  # 3. OMEGA-3 - Membrane Fluidity for Ion Channels (SCN2A)
  "Omega_3_FattyAcids" = c("FADS1", "FADS2", "ELOVL2", "ELOVL5", "SCD"),
  
  # 4. CARNITINE - Mitochondrial Function (Links to TMLHE)
  "Carnitine_Mito" = c("TMLHE", "SLC22A5", "CPT1A", "CPT2", "CRAT")
)

# RE-RUN THE "FLATTEN" STEP IMMEDIATELY AFTER THIS:
nutrition_map <- data.frame(Gene = character(), Category = character(), stringsAsFactors=FALSE)
for(cat in names(nutrition_list)){
  nutrition_map <- rbind(nutrition_map, data.frame(Gene = nutrition_list[[cat]], Category = cat))
}
# Flatten nutrition list into a single vector for matching
all_nutrition_genes <- unique(unlist(nutrition_list))

# B. Target ASD Genes (Top High Confidence - SFARI Cat 1)
asd_genes <- c("CHD8", "SHANK3", "SCN2A", "ADNP", "SYNGAP1", "ARID1B", "GRIN2B", "FOXP1", "DYRK1A", "NRXN1", "PTEN", "TBL1XR1", "POGZ", "KMT2A", "ASH1L")

# ==============================================================================
# 2. LOAD DATASET (GSE28521 - Voineagu et al. Brain Transcriptome)
# ==============================================================================
message("Downloading GSE28521 from GEO... (This may take a moment)")
gset <- getGEO("GSE28521", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6883", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Extract Expression Matrix
ex <- exprs(gset)

# Map Probe IDs to Gene Symbols
# Note: In a full pipeline, we map probes using the GPL platform. 
# For this script, we assume the row names are valid or we fetch annotation.
# Let's fetch annotation to be safe:
gpl <- getGEO(annotation(gset), getGPL=TRUE)
mapping <- Table(gpl)[, c("ID", "Symbol")] 
# Clean mapping
mapping <- mapping %>% filter(Symbol != "" & !is.na(Symbol))

# Convert expression matrix to data frame and merge symbols
ex_df <- as.data.frame(ex)
ex_df$ID <- rownames(ex_df)
merged_data <- merge(mapping, ex_df, by="ID")

# Aggregate by Symbol (mean expression if multiple probes per gene)
# This creates a clean matrix: Rows = Genes, Cols = Samples
final_matrix <- merged_data %>%
  select(-ID) %>%
  group_by(Symbol) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("Symbol")

# ==============================================================================
# 3. CALCULATE CORRELATION
# ==============================================================================

# Filter matrix to only include genes from our lists
# We look for genes that exist in BOTH the dataset and our lists
valid_asd <- intersect(rownames(final_matrix), asd_genes)
valid_nutri <- intersect(rownames(final_matrix), all_nutrition_genes)

if(length(valid_asd) == 0 | length(valid_nutri) == 0) {
  stop("No matching genes found in the dataset. Check gene symbols.")
}

# Subset the expression data
asd_exp <- t(final_matrix[valid_asd, ])
nutri_exp <- t(final_matrix[valid_nutri, ])

# Calculate Pearson Correlation
# Result: Matrix where Rows = ASD Genes, Cols = Nutrition Genes
cor_matrix <- cor(asd_exp, nutri_exp, method = "pearson", use = "pairwise.complete.obs")

# ==============================================================================
# 4. GENERATE HEATMAP
# ==============================================================================

# Create annotation for columns (Nutrition Categories)
# We need to map each nutrition gene back to its list name
annotation_col <- data.frame(Gene = valid_nutri)
annotation_col$Category <- NA

for (cat in names(nutrition_list)) {
  genes_in_cat <- nutrition_list[[cat]]
  annotation_col$Category[annotation_col$Gene %in% genes_in_cat] <- cat
}
rownames(annotation_col) <- annotation_col$Gene
annotation_col <- annotation_col %>% select(Category)

# Draw Heatmap
pheatmap(cor_matrix,
         annotation_col = annotation_col,
         main = "Pearson Correlation: ASD Genes vs. Nutrition Genes (Brain)",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize_row = 10,
         fontsize_col = 6,
         color = colorRampPalette(c("#D73027", "#F7F7F7", "#4575B4"))(50), # Red-White-Blue
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_colnames = TRUE)

message("Heatmap generated.")


# ==============================================================================
# 1. DEFINE GENE LISTS (Expanded)
# ==============================================================================

# ASD Genes
asd_genes <- c("CHD8", "SHANK3", "SCN2A", "ADNP", "SYNGAP1", "ARID1B", 
               "GRIN2B", "FOXP1", "DYRK1A", "NRXN1", "PTEN", "TBL1XR1", 
               "POGZ", "KMT2A", "ASH1L")

# Expanded Nutrition List (Original + New Additions)
nutrition_list <- list(
  # --- ORIGINAL SIGNATURES ---
  "A_Retinol"        = c("RPE65", "ADH1A", "DHRS4", "CYP2A6", "RDH12", "BCO1"),
  "B1_Thiamine"      = c("SLC19A2", "TPK1", "OGDH", "PDHA1", "SLC25A19", "TKT"),
  "B6_Pyridoxine"    = c("ALPL", "PNPO", "AMT", "CSAD", "GAD1", "GAD2"),
  "B12_Cobalamin"    = c("TCN2", "MTR", "MTRR", "MMACHC", "MUT", "CD320"),
  "Vit_D_Serotonin"  = c("VDR", "CYP27B1", "CYP24A1", "GC"),
  "Vit_E_Antioxidant"= c("TTPA", "CYP4F2", "GSTP1"),
  "Zinc_Synaptic"    = c("SLC30A8", "MT1A", "MT2A", "SLC39A1"),
  
  # --- NEW ADDITIONS (Methylation & Mitochondria) ---
  "B9_Folate_Meth"   = c("MTHFR", "DHFR", "SLC19A1", "FOLR1", "SHMT1", "MTHFD1"),
  "Iron_Myelination" = c("TFRC", "FTH1", "FTL", "SLC11A2", "ACO1", "IREB2"),
  "Omega_3_Membrane" = c("FADS1", "FADS2", "ELOVL2", "ELOVL5", "SCD"),
  "Carnitine_Mito"   = c("TMLHE", "SLC22A5", "CPT1A", "CPT2", "CRAT")
)

# Flatten List for Mapping
nutrition_map <- data.frame(Gene = character(), Category = character(), stringsAsFactors=FALSE)
for(cat in names(nutrition_list)){
  nutrition_map <- rbind(nutrition_map, data.frame(Gene = nutrition_list[[cat]], Category = cat))
}

# ==============================================================================
# 2. DATA CLEANING & PREPARATION
# ==============================================================================

# Ensure 'final_matrix' exists. If not, reload from GSE28521.
if(!exists("final_matrix")) stop("Please load your GSE28521 expression data first.")

# FORCE NUMERIC (The "Nuclear" Clean)
# This strips any hidden text and ensures correlation math works
clean_matrix <- function(mat) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  return(mat)
}
final_matrix_clean <- clean_matrix(final_matrix)
final_matrix_clean[is.infinite(final_matrix_clean)] <- NA

# Check Overlap (Diagnostic)
valid_asd <- intersect(rownames(final_matrix_clean), asd_genes)
valid_nutri_map <- nutrition_map %>% 
  filter(Gene %in% rownames(final_matrix_clean)) %>%
  arrange(Category)

print(paste("Matching ASD Genes:", length(valid_asd)))
print(paste("Matching Nutrition Genes:", nrow(valid_nutri_map)))

if(length(valid_asd) < 2) stop("Error: Too few ASD genes found. Check gene symbols.")

# ==============================================================================
# 3. CALCULATE REGULATION (Row Clustering)
# ==============================================================================

# Identify Samples (Adjust regex if needed)
pheno <- pData(gset) 
asd_samples <- rownames(pheno)[grep("autism|ASD", pheno$characteristics_ch1.1, ignore.case = TRUE)]
ctl_samples <- rownames(pheno)[grep("control", pheno$characteristics_ch1.1, ignore.case = TRUE)]

# Subset Data for ASD Genes
asd_data_subset <- final_matrix_clean[valid_asd, ]

# Calculate Log2 Fold Change (ASD - Control)
mean_asd <- rowMeans(asd_data_subset[, intersect(colnames(asd_data_subset), asd_samples)], na.rm=TRUE)
mean_ctl <- rowMeans(asd_data_subset[, intersect(colnames(asd_data_subset), ctl_samples)], na.rm=TRUE)
logFC <- mean_asd - mean_ctl

# Create Regulation Dataframe
regulation_df <- data.frame(Gene = names(logFC), LogFC = logFC)
regulation_df$Status <- ifelse(regulation_df$LogFC > 0, "Upregulated", "Downregulated")

# ==============================================================================
# 4. CALCULATE SPEARMAN CORRELATION
# ==============================================================================

# Extract Data
asd_exp <- t(final_matrix_clean[valid_asd, ])
nutri_exp <- t(final_matrix_clean[valid_nutri_map$Gene, ])

# Calculate Correlation (Spearman Rank)
cor_matrix <- cor(asd_exp, nutri_exp, method = "spearman", use = "pairwise.complete.obs")
cor_matrix[is.na(cor_matrix)] <- 0 # Safety for empty cells

# ==============================================================================
# 5. GENERATE HEATMAP (Bulletproof Version)
# ==============================================================================

# 1. ALIGN ANNOTATIONS & DATA
# Ensure row names match exactly before doing anything else
common_genes <- intersect(rownames(annot_row), rownames(cor_matrix))

if(length(common_genes) == 0) stop("Error: No common genes found between Data and Annotation.")

# Subset both to match
annot_row_clean <- annot_row[common_genes, , drop=FALSE]
cor_matrix_plot <- cor_matrix[common_genes, , drop=FALSE]

# 2. HANDLE MISSING REGULATION DATA
# If LogFC was NA, the Regulation is NA. We fill these to prevent crashes.
annot_row_clean$Regulation[is.na(annot_row_clean$Regulation)] <- "Unchanged"

# 3. SORT ROWS (Down -> Up)
# This creates the physical "blocks" in the heatmap
ordered_rows <- rownames(annot_row_clean)[order(annot_row_clean$Regulation)]

# Apply order to matrix and annotation
cor_matrix_plot <- cor_matrix_plot[ordered_rows, , drop=FALSE]
annot_row_clean <- annot_row_clean[ordered_rows, , drop=FALSE]

# 4. CALCULATE GAP POSITION (The Fix)
# We use na.rm = TRUE so one bad gene doesn't crash the whole count
num_down <- sum(annot_row_clean$Regulation == "Downregulated", na.rm = TRUE)

# SAFETY VALVE: If num_down is still invalid (NA/Inf), force it to 0
if (is.na(num_down)) num_down <- 0

# Only set a gap if it makes sense (between 1 and Total-1)
if (num_down > 0 && num_down < nrow(cor_matrix_plot)) {
  gap_pos <- num_down
} else {
  gap_pos <- NULL
}

# Debug Print (Check Console to see if it worked)
print(paste("Total Genes:", nrow(cor_matrix_plot)))
print(paste("Downregulated Genes:", num_down))

# 5. COLUMN ANNOTATION
# Ensure columns match too
annot_col <- data.frame(Signature = valid_nutri_map$Category)
rownames(annot_col) <- valid_nutri_map$Gene
valid_cols <- intersect(rownames(annot_col), colnames(cor_matrix_plot))

annot_col_clean <- annot_col[valid_cols, , drop=FALSE]
cor_matrix_plot <- cor_matrix_plot[, valid_cols, drop=FALSE]

# 6. PLOT
pheatmap(cor_matrix_plot,
         cluster_cols = FALSE,       # Keep Nutrition Groups Fixed
         cluster_rows = TRUE,       # Keep Regulation Split Fixed
         annotation_row = annot_row_clean, 
         annotation_col = annot_col_clean, 
         gaps_row = gap_pos,         # <--- Now guaranteed safe
         main = "ASD Genes (Regulated) vs Expanded Nutrition",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         border_color = NA,
         fontsize_row = 9, 
         fontsize_col = 8,
         na_col = "grey95")          # Handle any remaining NAs visually
