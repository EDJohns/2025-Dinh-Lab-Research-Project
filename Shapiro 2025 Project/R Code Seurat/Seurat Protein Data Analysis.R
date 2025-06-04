# Loading my Seurat Data

library(Seurat)
library(SeuratObject)
library(SeuratDisk)

# Load the Seurat object from the .rda file
load("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\R Data Objects\\HNSCC CosMx-selected\\CosMX_protein_annot_cluster.rda")

# Only run this if you're sure the object is named "CosMX_protein_annot_cluster"
SeuratObject <- get("CosMx_Protein_log_norm")

FeaturePlot(SeuratObject, features = c("CD3", "HER2"))

# loaded_names <- load("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\R Data Objects\\HNSCC CosMx-selected\\CosMX_protein_annot_cluster.rda")
# print(loaded_names)

