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
# RNA Counts (using k10_NhCoord20_anno instead of NhCoord20_anno)
# -----------------------------
meta_rna_k10 <- as.data.frame(table(CosMX_RNA_merged$k10_NhCoord20_anno, CosMX_RNA_merged$tma_fov_id))
colnames(meta_rna_k10) <- c("CellType", "FOV", "Freq")

meta_rna_k10 <- meta_rna_k10 %>%
  group_by(FOV) %>%
  mutate(Percentage = round(Freq / sum(Freq) * 100, 2)) %>%
  ungroup()

# -----------------------------
# Subset to DC_T_iCAF_Tumor_Mixed in RNA and Fibroblast_Immune_Mixed in Protein
# -----------------------------
rna_subset <- meta_rna_k10 %>%
  filter(CellType == "DC_T_iCAF_Tumor_Mixed")

protein_subset <- meta_protein %>%
  filter(CellType == "Fibroblast_Immune_Mixed")

# Merge RNA + Protein by FOV
rna_protein_subset <- merge(
  rna_subset[, c("FOV", "CellType", "Percentage")],
  protein_subset[, c("FOV", "CellType", "Percentage")],
  by = "FOV",
  suffixes = c("_RNA", "_Protein")
)


# -----------------------------
# RNA Counts (using k10_NhCoord20_anno instead of NhCoord20_anno)
# -----------------------------
meta_rna_k10 <- as.data.frame(table(CosMX_RNA_merged$k10_NhCoord20_anno, CosMX_RNA_merged$tma_fov_id))
colnames(meta_rna_k10) <- c("CellType", "FOV", "Freq")

meta_rna_k10 <- meta_rna_k10 %>%
  group_by(FOV) %>%
  mutate(Percentage = round(Freq / sum(Freq) * 100, 2)) %>%
  ungroup()

# -----------------------------
# Subset to CAF_Immune_Mixed or DC_T_iCAF_Tumor_Mixed in RNA
# and Fibroblast_Immune_Mixed in Protein
# -----------------------------
rna_subset <- meta_rna_k10 %>%
  filter(CellType %in% c("CAF_Immune_Mixed", "DC_T_iCAF_Tumor_Mixed"))

protein_subset <- meta_protein %>%
  filter(CellType == "Fibroblast_Immune_Mixed")

# Merge RNA + Protein by FOV
rna_protein_subset <- merge(
  rna_subset[, c("FOV", "CellType", "Percentage")],
  protein_subset[, c("FOV", "CellType", "Percentage")],
  by = "FOV",
  suffixes = c("_RNA", "_Protein")
)

# -----------------------------
# Plot: CAF_Immune_Mixed/DC_T_iCAF_Tumor_Mixed (RNA) vs Fibroblast_Immune_Mixed (Protein)
# -----------------------------
rna_protein_subset <- rna_protein_subset %>%
  mutate(TMA_Group = case_when(
    grepl("TMA_A", FOV) ~ "TMA_A",
    grepl("TMA_B", FOV) ~ "TMA_B",
    grepl("TMA_C", FOV) ~ "TMA_C",
    TRUE ~ "Other"
  ))

ggplot(rna_protein_subset, aes(x = Percentage_RNA, y = Percentage_Protein, color = TMA_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(color = TMA_Group), linetype = "dashed") +
  stat_poly_eq(
    aes(
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), 
      color = TMA_Group
    ),
    formula = y ~ x, 
    parse = TRUE, 
    show.legend = FALSE
  ) +
  theme_classic() +
  xlab("CAF_Immune_Mixed / DC_T_iCAF_Tumor_Mixed by FOV") +
  ylab("Fibroblast_Immune_Mixed by FOV in CosMx Protein") +
  ggtitle("CAF/DC_T_iCAF_Tumor_Mixed (RNA) vs Fibroblast_Immune_Mixed (Protein)") +
  scale_color_manual(values = c("TMA_A" = "#1b9e77", "TMA_B" = "#d95f02", "TMA_C" = "#7570b3")) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

########### PLOT 2: ALL IMMUNE CELLS #######################
# -----------------------------
# Define immune cell types
# -----------------------------
immune_types <- c("B_cells", "CD4+T_cells", "CD8+T_cells", 
                  "DCs", "Tregs", "Plasma_cells", 
                  "Monocytes", "Macrophages")

# -----------------------------
# Filter data to immune cells only
# -----------------------------
immune_data <- all_data %>%
  filter(CellType_RNA %in% immune_types) %>%
  mutate(TMA_Group = case_when(
    grepl("TMA_A", FOV_RNA) ~ "TMA_A",
    grepl("TMA_B", FOV_RNA) ~ "TMA_B",
    grepl("TMA_C", FOV_RNA) ~ "TMA_C",
    TRUE ~ "Other"
  ))

# -----------------------------
# Plot: All Immune Cells (RNA vs Protein)
# -----------------------------
ggplot(immune_data, aes(x = RNA_Percentage, y = Protein_Percentage, color = TMA_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(color = TMA_Group), linetype = "dashed") +
  stat_poly_eq(
    aes(
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"),
      color = TMA_Group
    ),
    formula = y ~ x,
    parse = TRUE,
    show.legend = FALSE
  ) +
  theme_classic() +
  xlab("Immune Cell Percentage by FOV in CosMx RNA data") +
  ylab("Immune Cell Percentage by FOV in CosMx Protein data") +
  ggtitle("Immune Cells Only - Cell Type Frequency per FOV (TMA A, B, C)") +
  scale_color_manual(values = c("TMA_A" = "#1b9e77",
                                "TMA_B" = "#d95f02",
                                "TMA_C" = "#7570b3")) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )



