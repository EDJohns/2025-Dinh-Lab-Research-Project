# Loading my Seurat Data

library(Seurat)
library(SeuratObject)
library(SeuratDisk)

# Load the Seurat object from the .rda file
load("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\R Data Objects\\HNSCC CosMx-selected\\CosMX_RNA_merged.rda")

loaded_names <- load("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\R Data Objects\\HNSCC CosMx-selected\\CosMX_RNA_merged.rda")
print(loaded_names)


mRNA_Seurat_Object <- get("CosMX_RNA_merged")

# Elbow plot of the PC, typically indicate different cell types
ElbowPlot(CosMX_RNA_merged)

# CosMX_RNA_merged <- FindNeighbors(CosMX_RNA_merged, dims = 1:20)
# CosMX_RNA_merged <- FindClusters(CosMX_RNA_merged, resolution = 1.2)

FeaturePlot(mRNA_Seurat_Object, features = c("CSF3R", "S100A8", "S100A9", "ELANE", "AZU1", "LCN2", "MPO"))

######################### Feature Plot of PMN CSF3R Gene ##############

FeaturePlot(mRNA_Seurat_Object, features = c("CSF3R"))


##################### Non Linear Dim Reduction (UMAP) #################
# Note that umap is 2 plot of cellular difference 

mRNA_Seurat_Object <- RunUMAP(mRNA_Seurat_Object, dims = 1:20)
DimPlot(mRNA_Seurat_Object, reduction = "umap")

##################### Plot Nuetorphil markers by PC Axes ########################


VlnPlot(mRNA_Seurat_Object, features = c("CSF3R"))
VlnPlot(mRNA_Seurat_Object, features = c("S100A8"))
VlnPlot(mRNA_Seurat_Object, features = c("S100A9"))

#################### Finding differently expressed Features ############

# find all markers of cluster 2
cluster2.markers <- FindMarkers(mRNA_Seurat_Object, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(mRNA_Seurat_Object, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
