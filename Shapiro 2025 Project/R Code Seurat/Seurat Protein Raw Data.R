# Seurat Raw Data Testing - This is to compare the raw data to the data that Athena provided to see if there is additional insights/changes I would make
# The data that Athena gave was already normalized. I want to compare that the raw data values


############################# Part 1: Loading the Data ############################# 

# Loading my Seurat Data

library(Seurat)
library(SeuratObject)

# Raw Protein Data Load - I am assuming there is no need to pull the 
seurat_object_raw <- readRDS("C:\\Users\\ejohns\\Downloads\\Protein_seurat_object.Rds")

# Load the Seurat object from the .rda file
load("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\R Data Objects\\HNSCC CosMx-selected\\CosMX_protein_annot_cluster.rda")
seurat_object_athena <- get("CosMx_Protein_log_norm")


############################# Part 2: Basic Exploration ############################# 

# Slot Names Examples
# Seurat Object (Toolbox)
# ├─ assays (Drawer: Expression data)
# ├─ meta.data (Drawer: Cell annotations)
# ├─ reductions (Drawer: PCA, UMAP, etc.)
# ├─ active.ident (Drawer: Cell clusters)
# ├─ graphs (Drawer: Network graphs)
slotNames(seurat_object_raw)
slotNames(seurat_object_athena)

# Find the dimensions (Features * Cells)
dim(seurat_object_raw)
dim(seurat_object_athena)

# List of available assays & expression values
Assays(seurat_object_raw)
Assays(seurat_object_athena)

# List of availble layers
Layers(seurat_object_raw[["negprobes"]])
Layers(seurat_object_raw[["RNA"]])
Layers(seurat_object_athena[["Protein"]])

# Peak at the expression values
seurat_object_raw[["negprobes"]]@data[1:5, 1:5]
seurat_object_athena[["Protein"]]@data[1:5, 1:5]

#Show Gene Names
rownames(seurat_object_athena[["Protein"]])
rownames(seurat_object_raw[["RNA"]])
rownames(seurat_object_raw[["negprobes"]])

#Showing expression data
GetAssayData(seurat_object_raw, assay = "RNA", layer = "counts")[1:15, 1:15]
GetAssayData(seurat_object_athena, assay = "Protein", layer = "counts")[1:15, 1:15]
