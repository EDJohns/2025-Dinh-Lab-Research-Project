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

# -----------------------------
# Plot 2: Immune Cell Types Only
# -----------------------------
immune_groups <- c("Immune_Mixed", "DC_T_B_NK_Centric", "Myeloid_Centric", "Fibroblast_Immune_Mixed")

immune_data <- all_data %>%
  filter(CellType_RNA %in% immune_groups)

ggplot(immune_data, aes(x = RNA_Percentage, y = Protein_Percentage, color = CellType_RNA)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, color = "black") +
  theme_classic() +
  xlab("Immune Cell Percentage by FOV mRNA (NhCoord20_anno)") +
  ylab("Immune Cell Percentage by FOV Protein (NhCoord20_anno)") +
  scale_color_manual(values = color_clusters) +
  labs(color = "Immune Cell Type") +
  ggtitle("Immune cell frequency per FOV (RNA vs Protein)") +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# -----------------------------
# Plot 3: Fibroblast_Immune_Mixed Only
# -----------------------------

#fibro_immune_data <- all_data %>%
#  filter(CellType_RNA == "Fibroblast_Immune_Mixed")

fibro_immune_data <- all_data %>%
  filter(CellType_RNA == "Fibroblast_Immune_Mixed" & grepl("TMA_C", Group))


ggplot(fibro_immune_data, aes(x = RNA_Percentage, y = Protein_Percentage)) +
  geom_point(color = "#E78AC3", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE, color = "black") +
  theme_classic() +
  xlab("Fibroblast_Immune_Mixed % by FOV mRNA") +
  ylab("Fibroblast_Immune_Mixed % by FOV Protein") +
  ggtitle("Fibroblast_Immune_Mixed frequency per FOV") +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# -----------------------------
# Plot 4: Fibroblast_Immune_Mixed split by TMA group
# -----------------------------

fibro_immune_data_tma <- all_data %>%
  filter(CellType_RNA == "Fibroblast_Immune_Mixed") %>%
  mutate(TMA_Group = case_when(
    grepl("TMA_A", Group) ~ "TMA_A",
    grepl("TMA_B", Group) ~ "TMA_B",
    grepl("TMA_C", Group) ~ "TMA_C",
    TRUE ~ "Other"
  ))

ggplot(fibro_immune_data_tma, aes(x = RNA_Percentage, y = Protein_Percentage, color = TMA_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(color = TMA_Group), linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"), color = TMA_Group), 
               formula = y ~ x, parse = TRUE, show.legend = FALSE) +
  theme_classic() +
  xlab("Fibroblast_Immune_Mixed by FOV in CosMx RNA") +
  ylab("Fibroblast_Immune_Mixed by FOV in CosMx Protein") +
  ggtitle("Fibroblast_Immune_Mixed frequency per FOV") +
  scale_color_manual(values = c("TMA_A" = "#1b9e77", "TMA_B" = "#d95f02", "TMA_C" = "#7570b3")) +
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
    Var1.x != "Tumor_cells",               # exclude tumor cells
    Var1.x %in% allowed_types              # only include immune cell types
  ) %>%
  mutate(TMA_Group = case_when(
    grepl("TMA_A", group) ~ "TMA_A",
    grepl("TMA_B", group) ~ "TMA_B",
    grepl("TMA_C", group) ~ "TMA_C",
    TRUE ~ "Other"
  ))

ggplot(filtered_data3, aes(x = Aggregated_Percentage, y = Percentage, color = TMA_Group)) +
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
