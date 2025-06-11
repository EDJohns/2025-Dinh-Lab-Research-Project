load("C:\\Users\\ejohns\\Documents\\GitHub\\2025-Dinh-Lab-Research-Project\\Shapiro 2025 Project\\Athena R Data Files\\CosMX_RNA_merged.rda")
load("C:\\Users\\ejohns\\Documents\\GitHub\\2025-Dinh-Lab-Research-Project\\Shapiro 2025 Project\\Athena R Data Files\\CosMX_protein_annot_cluster.rda")


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

library(knitr)
kable(head(meta_rna_percentages), caption = "Cell Type Percentages by Sample")
kable(head(meta_protein_percentages), caption = "Cell Type Percentages by Sample")
kable(tail(meta_rna_percentages), caption = "Cell Type Percentages by Sample")
kable(tail(meta_protein_percentages), caption = "Cell Type Percentages by Sample")

library(dplyr)

meta_rna_percentages %>%
  filter(Var2 == "TMA_A_FOV1") %>%
  arrange(desc(Percentage))  # Optional: sort by percentage descending

meta_rna_percentages %>%
  filter(Var2 == "TMA_A_FOV1") %>%
  kable(caption = "Cell Type Percentages in Sample TMA_A_FOV1")


meta_protein_percentages %>%
  filter(Var2 == "TMA_A_FOV1") %>%
  kable(caption = "Cell Type Percentages in Sample TMA_A_FOV1")

meta_rna_percentages %>%
  filter(Var2 == "TMA_C_FOV20") %>%
  kable(caption = "Cell Type Percentages in Sample TMA_C_FOV20")


meta_protein_percentages %>%
  filter(Var2 == "TMA_C_FOV20") %>%
  kable(caption = "Cell Type Percentages in Sample TMA_C_FOV20")


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






ggplot(all_data, aes(x=Aggregated_Percentage, y=Percentage, color = Var1.x)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, color = "black") +
  theme_classic() + 
  xlab('Cell Type Percentage by FOV in CosMx RNA data') + 
  ylab('Cell Type Percentage by FOV in CosMx Protein data') +
  scale_color_manual(values = color_clusters) +
  labs(color = "Cell Type") +
  ggtitle('Cell type frequency per FOV by dataset') +
  theme(
    axis.title = element_text(size = 16, face = "bold"),  # Larger axis labels
    axis.text = element_text(size = 14),  # Larger axis text
    legend.title = element_text(size = 14, face = "bold"),  # Larger legend title
    legend.text = element_text(size = 12),  # Larger legend text
    legend.key.size = unit(1, "cm"),  # Larger legend keys
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Larger, centered title
  )

ggsave('cosmx_celltypecomparison_scatterplot.pdf', dpi = 900, width = 9, height = 6)

kable(all_data, caption = "Cell Type Percentages by Sample")
kable(head(meta_protein_percentages), caption = "Cell Type Percentages by Sample")


CosMX_RNA_merged <- CosMX_RNA_merged[,grep('Tumor', CosMX_RNA_merged$annot_cluster2, invert = T)]
CosMx_Protein_log_norm <- CosMx_Protein_log_norm[,grep('Tumor', CosMx_Protein_log_norm$merged_annot_cluster, invert = T)]

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






ggplot(all_data, aes(x=Aggregated_Percentage, y=Percentage, color = Var1.x)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, color = "black") +
  theme_classic() + 
  xlab('Cell Type Percentage by FOV in CosMx RNA data') + 
  ylab('Cell Type Percentage by FOV in CosMx Protein data') +
  scale_color_manual(values = color_clusters) +
  labs(color = "Cell Type") +
  ggtitle('Cell type frequency per FOV by dataset') +
  theme(
    axis.title = element_text(size = 16, face = "bold"),  # Larger axis labels
    axis.text = element_text(size = 14),  # Larger axis text
    legend.title = element_text(size = 14, face = "bold"),  # Larger legend title
    legend.text = element_text(size = 12),  # Larger legend text
    legend.key.size = unit(1, "cm"),  # Larger legend keys
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Larger, centered title
  )

ggsave('NOTUMORS_cosmx_celltypecomparison_scatterplot.pdf', dpi = 900, width = 9, height = 6)

kable(all_data, caption = "Cell Type Percentages by Sample")

