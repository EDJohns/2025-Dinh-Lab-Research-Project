# -----------------------------
# Load Data
# -----------------------------
load("C:\\Users\\evanj\\OneDrive\\Documents\\Shapiro Data Files\\Rstudio RDA Files\\CosMX_RNA_merged.rda")
load("C:\\Users\\evanj\\OneDrive\\Documents\\Shapiro Data Files\\Rstudio RDA Files\\CosMX_protein_annot_cluster.rda")

# -----------------------------
# Packages
# -----------------------------
packages <- c("dplyr", "tidyverse", "ggpmisc", "stringr")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}
library(dplyr)
library(tidyverse)
library(ggpmisc)
library(stringr)

# -----------------------------
# Color Palette
# -----------------------------
color_clusters <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", 
  "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", 
  "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", 
  "#666666", "#999999", "#AA8282", "#D4B7B7", "#8600BF", "#BA5CE3", 
  "#808000","#AEAE5C", "#1E90FF", "#00BFFF", "#56FF0D", "#FFFF00"
)

# -----------------------------
# Mapping function (RNA → Protein)
# -----------------------------
map_rna_to_protein <- function(x) {
  case_when(
    x == "Endothelial_Centric"     ~ "Endothelial_Centric",
    x == "Bcell_Plasma_Centric"    ~ "DC_T_B_NK_Centric",
    x == "Immune_Mixed_Mac_Neg"    ~ "Immune_Mixed",
    x == "Fibroblast_Centric"      ~ "Fibroblast_Centric",
    x == "Fibroblast_Immune_Mixed" ~ "Fibroblast_Immune_Mixed",
    x == "Macrophage_Centric"      ~ "Myeloid_Centric",
    x == "Immune_Mixed_CD4Thi"     ~ "Immune_Mixed",
    x == "CD8T_Tumor_Mixed"        ~ "Immune_Mixed",
    x == "Tumor_Epi_Centric"       ~ "Tumor_Epi_Centric",
    TRUE ~ x
  )
}

# -----------------------------
# RNA Counts → Percentages
# -----------------------------
meta_rna <- as.data.frame(table(CosMX_RNA_merged$NhCoord20_anno, CosMX_RNA_merged$tma_fov_id))
colnames(meta_rna) <- c("CellType", "FOV", "Freq")

meta_rna <- meta_rna %>%
  group_by(FOV) %>%
  mutate(Percentage = round(Freq / sum(Freq) * 100, 2)) %>%
  ungroup() %>%
  mutate(CellType_mapped = map_rna_to_protein(CellType))

# -----------------------------
# Protein Counts → Percentages
# -----------------------------
meta_protein <- as.data.frame(table(CosMx_Protein_log_norm$NhCoord20_anno, CosMx_Protein_log_norm$tma_fov_id))
colnames(meta_protein) <- c("CellType", "FOV", "Freq")

meta_protein <- meta_protein %>%
  group_by(FOV) %>%
  mutate(Percentage = round(Freq / sum(Freq) * 100, 2)) %>%
  ungroup()

# -----------------------------
# Merge RNA + Protein by CellType + FOV
# -----------------------------
meta_rna$Group <- paste(meta_rna$CellType_mapped, meta_rna$FOV, sep=":")
meta_protein$Group <- paste(meta_protein$CellType, meta_protein$FOV, sep=":")

all_data <- merge(
  meta_rna[, c("Group", "CellType_mapped", "FOV", "Percentage")],
  meta_protein[, c("Group", "CellType", "FOV", "Percentage")],
  by = "Group"
)

colnames(all_data) <- c("Group", "CellType_RNA", "FOV_RNA", "RNA_Percentage", 
                        "CellType_Protein", "FOV_Protein", "Protein_Percentage")

# -----------------------------
# Plot: RNA vs Protein Frequencies
# -----------------------------
ggplot(all_data, aes(x = RNA_Percentage, y = Protein_Percentage, color = CellType_RNA)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, color = "black") +
  theme_classic() +
  xlab("Cell Type Percentage by FOV mRNA (NhCoord20_anno)") +
  ylab("Cell Type Percentage by FOV Protein (NhCoord20_anno)") +
  scale_color_manual(values = color_clusters) +
  labs(color = "Cell Type") +
  ggtitle("Cell type frequency per FOV (RNA vs Protein)") +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

######################### Immune Cells Only (TMA A, B, C together) #################

allowed_types <- c("B_cells", "CD4+T_cells", "CD8+T_cells", 
                   "DCs", "Tregs", "Plasma_cells", 
                   "Monocytes", "Macrophages")

filtered_data3 <- all_data %>%
  filter(
    CellType_RNA %in% allowed_types         # only include immune cell types
  ) %>%
  mutate(TMA_Group = case_when(
    grepl("TMA_A", Group) ~ "TMA_A",
    grepl("TMA_B", Group) ~ "TMA_B",
    grepl("TMA_C", Group) ~ "TMA_C",
    TRUE ~ "Other"
  ))

ggplot(filtered_data3, aes(x = RNA_Percentage, y = Protein_Percentage, color = TMA_Group)) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_smooth(method = "lm", se = FALSE, aes(color = TMA_Group), linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"), color = TMA_Group), 
               formula = y ~ x, parse = TRUE, show.legend = FALSE) +
  theme_classic() + 
  xlab('Immune Cell Percentage by FOV in CosMx RNA data') + 
  ylab('Immune Cell Percentage by FOV in CosMx Protein data') +
  scale_color_manual(values = c("TMA_A" = "#1b9e77", "TMA_B" = "#d95f02", "TMA_C" = "#7570b3")) +
  labs(color = "TMA Group") +
  ggtitle('Immune Cells Only - Cell type frequency per FOV (TMA A, B, C)') +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

library(Seurat)
# -----------------------------
# Option 1: Sample-level averages
# -----------------------------
# Genes of interest
genes_of_interest <- c("IFNG","CXCR3","CXCR6", "CXCL9", "CXCL10", "CXCL16")

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
# Option 1: Sample-level averages (Response == "R")
# -----------------------------

# Genes of interest
genes_of_interest <- c("IFNG","CXCR3","CXCR6", "CXCL9", "CXCL10", "CXCL16")

# Extract expression matrix for those genes + metadata
expr_mat <- FetchData(CosMX_RNA_merged, vars = c(genes_of_interest, "tma_fov_id", "Response"))

# Keep only samples where Response == "R"
expr_mat <- expr_mat %>% filter(Response == "NR")

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
  ggtitle("IFNG Correlation (Cell Type, Response = R)")


# -----------------------------
# Generate 12 heatmaps (6 cell types x 2 responses)
# -----------------------------

library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)

# Genes of interest
genes_of_interest <- c("IFNG","CXCR3","CXCR6", "CXCL9", "CXCL10", "CXCL16")

# Extract expression + metadata
expr_mat <- FetchData(
  CosMX_RNA_merged, 
  vars = c(genes_of_interest, "tma_fov_id", "Response", "cell_type")
)

# Define cell types and responses
cell_types <- c("DCs", "Monocytes", "Macrophages", "Tcells_CD4+", "Tcells_CD8+", "Tregs")
responses  <- c("R", "NR")

# Function to create one heatmap
make_heatmap <- function(df, cell_type, resp) {
  df_sub <- df %>%
    filter(Response == resp, cell_type == !!cell_type) %>%
    group_by(tma_fov_id) %>%
    summarise(across(all_of(genes_of_interest), mean, na.rm = TRUE), .groups = "drop")
  
  # If not enough data, skip
  if (nrow(df_sub) < 2) return(NULL)
  
  cor_matrix <- cor(df_sub[, genes_of_interest], use = "pairwise.complete.obs", method = "pearson")
  cor_df <- melt(cor_matrix)
  
  ggplot(cor_df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    ggtitle(paste0("IFNG Correlation (", cell_type, ", Response = ", resp, ")"))
}

# Generate list of plots
plots <- map(cell_types, function(ct) {
  map(responses, function(resp) make_heatmap(expr_mat, ct, resp))
}) %>% flatten()

# Example: view first plot
plots[[1]]

# If you want to arrange them all in a grid:
library(patchwork)
wrap_plots(plots, ncol = 3)


# -----------------------------
# One multidimensional figure (faceted heatmap grid)
# -----------------------------

library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)

# Genes of interest
genes_of_interest <- c("IFNG","CXCR3","CXCR6", "CXCL9", "CXCL10", "CXCL16")

# Extract expression + metadata
expr_mat <- FetchData(
  CosMX_RNA_merged, 
  vars = c(genes_of_interest, "tma_fov_id", "Response", "cell_type")
)

# Define cell types and responses
cell_types <- c("DCs", "Monocytes", "Macrophages", "Tcells_CD4+", "Tcells_CD8+", "Tregs")
responses  <- c("R", "NR")

# Function to compute correlation data for one subset
compute_cor_df <- function(df, cell_type, resp) {
  df_sub <- df %>%
    filter(Response == resp, cell_type == !!cell_type) %>%
    group_by(tma_fov_id) %>%
    summarise(across(all_of(genes_of_interest), mean, na.rm = TRUE), .groups = "drop")
  
  if (nrow(df_sub) < 2) return(NULL)
  
  cor_matrix <- cor(df_sub[, genes_of_interest], use = "pairwise.complete.obs", method = "pearson")
  cor_df <- melt(cor_matrix)
  cor_df$cell_type <- cell_type
  cor_df$Response <- resp
  cor_df
}

# Build full correlation dataframe
cor_all <- map_dfr(cell_types, function(ct) {
  map_dfr(responses, function(resp) compute_cor_df(expr_mat, ct, resp))
})

# Faceted heatmap: rows = cell types, cols = response (R/NR)
ggplot(cor_all, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 2, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  facet_grid(cell_type ~ Response) +
  ggtitle("IFNG-related correlations across cell types and responses") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )

