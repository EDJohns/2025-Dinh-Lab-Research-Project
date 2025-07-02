# Loading my Seurat Data

library(Seurat)
library(SeuratObject)
library(SeuratDisk)

# Try Loading in the data
seurat_object <- readRDS("C:\\Users\\ejohns\\Downloads\\seurat_object.Rds")


# Load the Seurat object from the .rda file
load("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\R Data Objects\\HNSCC CosMx-selected\\CosMX_protein_annot_cluster.rda")

"C:\\Users\\ejohns\\Downloads\\Protein_seurat_object.Rds"

loaded_names <- load("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\R Data Objects\\HNSCC CosMx-selected\\CosMX_protein_annot_cluster.rda")
print(loaded_names)

# Only run this if you're sure the object is named "CosMX_protein_annot_cluster"
SeuratObject <- get("CosMx_Protein_log_norm")

FeaturePlot(SeuratObject, features = c("CD3", "HER2"))

ElbowPlot(SeuratObject)

a#freature plot
FeaturePlot(SeuratObject, features = c("CD15"))
