load("C:\\Users\\evanj\\OneDrive\\Documents\\Shapiro Data Files\\Rstudio RDA Files\\CosMX_RNA_merged.rda")

load("C:\\Users\\evanj\\OneDrive\\Documents\\Shapiro Data Files\\Rstudio RDA Files\\CosMX_protein_annot_cluster.rda")


# Install packages if not already installed
packages <- c("dplyr", "tidyverse")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}


# Load the packages
library(dplyr)
library(tidyverse)
library(ggpmisc)


color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#FF7F00", "#FDB462", 
                    "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", 
                    "#7570B3", "#BEAED4", "#666666", "#999999", "#AA8282", "#D4B7B7", "#8600BF", "#BA5CE3", 
                    "#808000","#AEAE5C", "#1E90FF", "#00BFFF", "#56FF0D", "#FFFF00")

meta_rna <- as.data.frame(table(CosMX_RNA_merged$annot_cluster2, CosMX_RNA_merged$tma_fov_id))
meta_protein <- as.data.frame(table(CosMx_Protein_log_norm$merged_annot_cluster, CosMx_Protein_log_norm$tma_fov_id))

# calculate percentages by FOV for the RNA FOVs
# Step 2: Group by sample (Var2) and calculate percentages
meta_rna_percentages <- meta_rna %>%
  group_by(Var2) %>%
  mutate(Percentage = Freq / sum(Freq) * 100) %>%
  ungroup()
# Step 3: Round the percentages to two decimal places (optional)
meta_rna_percentages$Percentage <- round(meta_rna_percentages$Percentage, 2)
# Step 4: View the first few rows of the new dataframe
head(meta_rna_percentages)
# Sanity check
sanity_check <- meta_rna_percentages %>%
  group_by(Var2) %>%
  summarize(
    Total_Percentage = sum(Percentage),
    Is_100 = abs(Total_Percentage - 100) < 0.05
  ) %>%
  ungroup()
# View the results
print(sanity_check)
# Check if all samples sum to approximately 100%
all_correct <- all(sanity_check$Is_100)
cat("All samples sum to 100% (within rounding error):", all_correct, "\n")
# If not all correct, view the problematic samples
if (!all_correct) {
  cat("Samples not summing to 100% (within rounding error):\n")
  print(sanity_check[!sanity_check$Is_100, ])
}






# calculate percentages by FOV for the protein FOVs
# Step 2: Group by sample (Var2) and calculate percentages
meta_protein_percentages <- meta_protein %>%
  group_by(Var2) %>%
  mutate(Percentage = Freq / sum(Freq) * 100) %>%
  ungroup()
# Step 3: Round the percentages to two decimal places (optional)
meta_protein_percentages$Percentage <- round(meta_protein_percentages$Percentage, 2)
# Step 4: View the first few rows of the new dataframe
head(meta_protein_percentages)
# Sanity check
sanity_check <- meta_protein_percentages %>%
  group_by(Var2) %>%
  summarize(
    Total_Percentage = sum(Percentage),
    Is_100 = abs(Total_Percentage - 100) < 0.05
  ) %>%
  ungroup()
# View the results
print(sanity_check)
# Check if all samples sum to approximately 100%
all_correct <- all(sanity_check$Is_100)
cat("All samples sum to 100% (within rounding error):", all_correct, "\n")
# If not all correct, view the problematic samples
if (!all_correct) {
  cat("Samples not summing to 100% (within rounding error):\n")
  print(sanity_check[!sanity_check$Is_100, ])
}




unique(meta_rna_percentages$Var1)
# [1] Bcells                                 DCs                                    Endothelials                           Fibroblast_apCAFs                     
# [5] Fibroblast_myCAFs                      Fibroblast_VEGFA_inflammatorylike_CAFs Macrophages                            Monocytes                             
# [9] Plasma                                 Tcells_CD4+                            TcellsCD8+                             Tregs                                 
# [13] Tumors_HIF1A+                          Tumors_KRT17+                          Tumors_KRT19+                          Tumors_Others 
unique(meta_protein_percentages$Var1)
# [1] B_cells           CD4+T_cells       CD8+T_cells       DCs               Endothelial_cells Fibroblasts/SMCs  Macrophages       Monocytes         Neutrophils       NK_cells         
# [11] Plasma_cells      Tregs             Tumor_cells     





library(stringr)

# Create a new column for the broader categories
meta_rna_aggregated <- meta_rna_percentages %>%
  mutate(
    Broad_Category = case_when(
      str_starts(Var1, "Tumors") ~ "Tumor_cells",
      str_starts(Var1, "Fibroblast") ~ "Fibroblasts/SMCs",
      TRUE ~ as.character(Var1)  # Keep other categories as they are
    )
  )
# Aggregate the percentages for the broader categories
meta_rna_aggregated <- meta_rna_aggregated %>%
  group_by(Var2, Broad_Category) %>%
  summarize(
    Aggregated_Percentage = sum(Percentage),
    .groups = "drop"
  ) %>%
  rename(Var1 = Broad_Category)
# Round the aggregated percentages to two decimal places
meta_rna_aggregated$Aggregated_Percentage <- round(meta_rna_aggregated$Aggregated_Percentage, 2)
# View the first few rows of the new dataframe
head(meta_rna_aggregated)
# Perform a sanity check on the aggregated data
sanity_check <- meta_rna_aggregated %>%
  group_by(Var2) %>%
  summarize(
    Total_Percentage = sum(Aggregated_Percentage),
    Is_100 = abs(Total_Percentage - 100) < 0.05
  ) %>%
  ungroup()
print(sanity_check)
all_correct <- all(sanity_check$Is_100)
cat("All samples sum to 100% (within rounding error):", all_correct, "\n")
if (!all_correct) {
  cat("Samples not summing to 100% (within rounding error):\n")
  print(sanity_check[!sanity_check$Is_100, ])
}


meta_rna_aggregated <- meta_rna_aggregated %>%
  mutate(Var1 = case_when(
    Var1 == "Bcells" ~ "B_cells",
    Var1 == "Endothelials" ~ "Endothelial_cells",
    Var1 == "Plasma" ~ "Plasma_cells",
    Var1 == "Tcells_CD4+" ~ "CD4+T_cells",
    Var1 == "TcellsCD8+" ~ "CD8+T_cells",
    TRUE ~ as.character(Var1)  # Keep other categories as they are
  ))
# View the first few rows to check the changes
head(meta_rna_aggregated)
# Optional: Check unique values in Var1 to ensure all changes were made
unique(meta_rna_aggregated$Var1)
# look at the overlap between the dataset labels
intersect(meta_rna_aggregated$Var1, meta_protein_percentages$Var1)
# remove cell types not present in both groups so it is equal
meta_protein_percentages_mini <- meta_protein_percentages[grep('NK_cells|Neutrophils', meta_protein_percentages$Var1, invert = T), c('Var1', 'Var2', 'Percentage')]
# making a "group" column to merge by
meta_rna_aggregated$group <- paste(meta_rna_aggregated$Var1, meta_rna_aggregated$Var2, sep = ':')
meta_protein_percentages_mini$group <- paste(meta_protein_percentages_mini$Var1, meta_protein_percentages_mini$Var2, sep = ':')
# merge both RNa and protein data together
all_data <- merge(meta_rna_aggregated, meta_protein_percentages_mini, by = 'group')


# -----------------------------
# Sample Average Correlations
# -----------------------------

library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
# Genes of interest
genes_of_interest <- c(
  # Core signaling
  "IFNG", "IFNGR1", "IFNGR2",
  "JAK1", "JAK2",
  "STAT1", "STAT3", "STAT4",
  
  # Antigen processing & presentation
  "B2M", "TAP1", "TAP2", "CIITA",
  "HLA-A", "HLA-B", "HLA-C", "HLA-E",
  "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
  "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
  
  # IFN-γ–induced chemokines & cytokines
  "CXCL9", "CXCL10", "CXCL16",
  "CCL2", "CCL5", "CCL7", "CCL8",
  
  # Regulators / effectors
  "IRF3", "IRF4", "ICAM1"
)
# Extract expression matrix for those genes
expr_mat <- FetchData(CosMX_RNA_merged, vars = genes_of_interest)

# Add metadata (FOV IDs)
expr_mat$tma_fov_id <- CosMX_RNA_merged@meta.data$tma_fov_id

# Compute average expression per FOV
rna_avg_fov <- expr_mat %>%
  group_by(tma_fov_id) %>%
  summarise(across(all_of(genes_of_interest), mean, na.rm = TRUE))

# Correlation matrix
cor_matrix <- cor(rna_avg_fov[, genes_of_interest], use = "pairwise.complete.obs", method = "pearson")

# Plot heatmap
library(reshape2)
cor_df <- melt(cor_matrix)

ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  ggtitle("Correlation of CXCR3, CXCL9, CXCL16 vs IFNG (Sample-level averages)")

# -----------------------------
# Monocyte & Macrophage correlations
# -----------------------------

library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)

# Genes and cell types of interest
genes_of_interest <- c("CXCR3", "CXCL9", "CXCL16", "IFNG")
celltypes_of_interest <- c("Monocytes")  # adjust if names differ in your metadata

# Pull expression values for these genes
expr_mat_ct <- FetchData(CosMX_RNA_merged, vars = genes_of_interest)

# Add metadata (FOV + CellType from your object)
expr_mat_ct$tma_fov_id <- CosMX_RNA_merged@meta.data$tma_fov_id
expr_mat_ct$CellType   <- CosMX_RNA_merged@meta.data$annot_cluster2

# Restrict to monocytes + macrophages
expr_mat_ct <- expr_mat_ct %>%
  filter(CellType %in% celltypes_of_interest)

# Compute average expression per gene per FOV *within each cell type*
rna_avg_ct <- expr_mat_ct %>%
  group_by(CellType, tma_fov_id) %>%
  summarise(across(all_of(genes_of_interest), mean, na.rm = TRUE), .groups = "drop")

# Loop through cell types to generate correlation heatmaps
for(ct in celltypes_of_interest) {
  ct_data <- rna_avg_ct %>% filter(CellType == ct)
  
  cor_matrix_ct <- cor(ct_data[, genes_of_interest],
                       use = "pairwise.complete.obs",
                       method = "pearson")
  
  cor_df_ct <- melt(cor_matrix_ct)
  
  p <- ggplot(cor_df_ct, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    ggtitle(paste("Correlation of CXCR3, CXCL9, CXCL16 vs IFNG (", ct, " averages)", sep = ""))
  
  print(p)
}
# -----------------------------
# -----------------------------
# Correlation of genes vs IFNG
# -----------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# Genes of interest (exclude IFNG itself)
genes_of_interest <- c(
  "IFNGR1", "IFNGR2",
  "JAK1", "JAK2",
  "STAT1", "STAT3", "STAT4",
  "B2M", "TAP1", "TAP2", "CIITA",
  "HLA-A", "HLA-B", "HLA-C", "HLA-E",
  "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
  "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
  "CXCL9", "CXCL10", "CXCL16",
  "CCL2", "CCL5", "CCL7", "CCL8",
  "IRF3", "IRF4", "ICAM1"
)

# Extract expression matrix for those genes + IFNG
expr_mat <- FetchData(CosMX_RNA_merged, vars = c("IFNG", genes_of_interest))

# Add metadata (FOV IDs)
expr_mat$tma_fov_id <- CosMX_RNA_merged@meta.data$tma_fov_id

# Compute average expression per FOV
rna_avg_fov <- expr_mat %>%
  group_by(tma_fov_id) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Compute correlations vs IFNG
cor_values <- sapply(genes_of_interest, function(gene) {
  cor(rna_avg_fov[[gene]], rna_avg_fov$IFNG, use = "pairwise.complete.obs", method = "pearson")
})

cor_df <- data.frame(
  Gene = names(cor_values),
  Correlation = cor_values
) %>%
  arrange(desc(Correlation)) %>%
  mutate(Gene = factor(Gene, levels = Gene))  # preserve sorted order

# Plot correlations as bar/line
ggplot(cor_df, aes(x = Gene, y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  xlab("Gene") +
  ylab("Pearson Correlation with IFNG") +
  ggtitle("Correlation of Genes vs IFNG (all cell types)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# -----------------------------
# Correlation of genes vs IFNG (Monocytes & Macrophages)
# -----------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)

# Genes of interest (exclude IFNG itself)
genes_of_interest <- c(
  "IFNGR1", "IFNGR2",
  "JAK1", "JAK2",
  "STAT1", "STAT3", "STAT4",
  "B2M", "TAP1", "TAP2", "CIITA",
  "HLA-A", "HLA-B", "HLA-C", "HLA-E",
  "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
  "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
  "CXCL9", "CXCL10", "CXCL16",
  "CCL2", "CCL5", "CCL7", "CCL8",
  "IRF3", "IRF4", "ICAM1"
)

# Cell types of interest
celltypes_of_interest <- c("Macrophages")  # adjust if needed

# Pull expression for genes + metadata
expr_mat <- FetchData(CosMX_RNA_merged, vars = c("IFNG", genes_of_interest))
expr_mat$CellType <- CosMX_RNA_merged@meta.data$annot_cluster2
expr_mat$tma_fov_id <- CosMX_RNA_merged@meta.data$tma_fov_id

# Restrict to monocytes & macrophages
expr_mat_ct <- expr_mat %>%
  filter(CellType %in% celltypes_of_interest)

# Compute average expression per FOV for each cell type
rna_avg_ct <- expr_mat_ct %>%
  group_by(CellType, tma_fov_id) %>%
  summarise(across(c("IFNG", all_of(genes_of_interest)), mean, na.rm = TRUE), .groups = "drop")

# Compute correlations vs IFNG per cell type
cor_df_list <- lapply(celltypes_of_interest, function(ct) {
  ct_data <- rna_avg_ct %>% filter(CellType == ct)
  
  cor_values <- sapply(genes_of_interest, function(gene) {
    cor(ct_data[[gene]], ct_data$IFNG, use = "pairwise.complete.obs", method = "pearson")
  })
  
  data.frame(
    Gene = names(cor_values),
    Correlation = cor_values,
    CellType = ct
  )
})

cor_df <- do.call(rbind, cor_df_list) %>%
  group_by(CellType) %>%
  arrange(desc(Correlation), .by_group = TRUE) %>%
  mutate(Gene = factor(Gene, levels = Gene))  # preserve sorted order

# Plot correlations per cell type
ggplot(cor_df, aes(x = Gene, y = Correlation, fill = CellType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  xlab("Gene") +
  ylab("Pearson correlation with IFNG") +
  ggtitle("Correlation of genes vs IFNG (Macrophages)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Monocytes" = "#1b9e77", "Macrophages" = "#d95f02"))
