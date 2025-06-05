library(dplyr)
library(Seurat)
library(patchwork)

#Create a PBMC object
#pbmc.data <- CreateSeuratObject()

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/ejohns/Documents/Shapiro Data Files/Seurat Tutorial/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Print out of the first 30 cells. The C selects particular rows, 
# var.[c(rows desired)/ range(1:X), c(columns desired) range(1:X)]
pbmc.data[c("CD3D","TCL1A", "MS4A1"), 1:30]
#Note: The dots on the printout present zeros makes non zero numbers more clear

#################### Further Analysis ######################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 15)

#creating a violin plot
# VlnPlot Summary:
# VlnPlot() generates violin plots to visualize the distribution of a specific gene(s) 
# or metadata feature across clusters or groups of cells.
# Commonly used to inspect gene expression (e.g., QC metrics like nFeature_RNA, percent.mt).
# 
# Key arguments:
# - object: Seurat object
# - features: character vector of features to plot
# - group.by: grouping variable (e.g., clusters, sample ID)
# - pt.size: size of individual points (set to 0 to hide)
# - assay: which assay to use (e.g., "RNA")
#
# Example:
# VlnPlot(seurat_obj, features = "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0.1)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)


# FeatureScatter Summary:
# -----------------------
# FeatureScatter() creates a scatter plot to visualize the relationship between two 
# features (typically gene expression or QC metrics) across all cells.
# Commonly used for QC (e.g., nCount_RNA vs. percent.mt) to identify outliers.

# Key arguments:
# - object: Seurat object
# - feature1: x-axis feature (e.g., "nCount_RNA")
# - feature2: y-axis feature (e.g., "percent.mt")
# - cols: colors for points
# - pt.size: size of points in the scatter plot

# Example:
# FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot3 <- FeatureScatter(pbmc, feature1="nCount_RNA",feature2="nFeature_RNA")
plot4 <- FeatureScatter(pbmc, feature1="nCount_RNA", feature2 = "percent.mt")
plot3 + plot4 #Trans

# Creating subset for analysis
pbmc_subset <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#New VinPlot of the new and old data sets
plotVin1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
plotVin2 <- VlnPlot(pbmc_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
plotVin1
plotVin2

############################### Data Normalization #############################

# Normilze data function turns it into a log function and store the infomration in pbmc[["RNA"]]$data.
## | Part            | Meaning                                                            |
##  | --------------- | ------------------------------------------------------------------ |
##   | `pbmc`          | A Seurat object                                                    |
##   | `pbmc[["RNA"]]` | Gets the `"RNA"` Assay object from the list of assays              |
##  | `$data`         | Gets the `data` slot (log-normalized matrix) from the Assay object |

pbmc_norm <-NormalizeData(pbmc_subset, normalization.method = "LogNormalize",scale.factor=10000)

# plotVin1 <- VlnPlot(pbmc_norm[["RNA"]]$data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
# plotVin2 <- VlnPlot(pbmc_norm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
# plotVin1+plotVin2


########################### Highly Variable Features Analysis ##################

# Program to for RNA expression that is high in some cells but otherwise low in others. 
# I immagine this could be useful in identifying cells

pbmc_norm <- FindVariableFeatures(pbmc_norm, selection.method = "vst", nfeatures = 2000) #I believe 200 is the ubme of genes analyzed

# Print out top 12 most variable genes cell to cell
top12 <- head(VariableFeatures(pbmc_norm),20)

# plot variable features with and without labels. Note that both plots are custom to seurat
plot1 <- VariableFeaturePlot(pbmc_norm)
plot2 <- LabelPoints(plot=plot1,points=top12,repel=TRUE)
plot1+plot2

###################### Scaling the Data ############################

# shifts expression so mean is 0 for each gene. Changes variation is 1.
# Also stored pbmc[["RNA"]]$scale.data
all.genes <- rownames(pbmc_norm)
pbmc_norm <- ScaleData(pbmc_norm, features = all.genes)
pbmc_norm <- ScaleData(pbmc_norm, vars.to.regress = "percent.mt")

################## Linear Dimensional Reduction ####################

# Basically creates axes of variation. 
  ## | **PC1** | Direction of **greatest variance** in the data |
  ## | **PC2** | Second-most variance, orthogonal to PC1        |
  ## | **PC3** | Third-most variance, orthogonal to PC1 and PC2 |
  
pbmc_norm <- RunPCA(pbmc_norm, features = VariableFeatures(object = pbmc_norm))
# Examine and visualize PCA results a few different ways
print(pbmc_norm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc_norm, dims = 1:2, reduction = "pca")
DimPlot(pbmc_norm, reduction = "pca") + NoLegend()

DimHeatmap(pbmc_norm, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc_norm, dims = 1:15, cells = 500, balanced = TRUE)

######################## Dataset Dimensionality #####################

# Allows you to see how predictivie each PC is. In this case top 8 have most of the varience
ElbowPlot(pbmc_norm)

####################### Cluster The Cells ###########################

pbmc_norm <- FindNeighbors(pbmc_norm, dism = 1:20)
pbmc_norm <- FindClusters(pbmc_norm, resolution = 0.5)
head(Idents(pbmc_norm),20)

##################### Non Linear Dim Reduction (UMAP) #################
# Note that umap is 2 plot of cellular difference 

pbmc_norm <- RunUMAP(pbmc_norm, dims = 1:20)
DimPlot(pbmc_norm, reduction = "umap")

#################### Finding differently expressed Features ############

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc_norm, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc_norm, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# 
cluster0.markers <- FindMarkers(pbmc_norm, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc_norm, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc_norm, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#freature plot
FeaturePlot(pbmc_norm, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

## Assigning Cell IDs to different data clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc_norm)
pbmc_norm <- RenameIdents(pbmc_norm, new.cluster.ids)
DimPlot(pbmc_norm, reduction = "umap", label = TRUE, pt.size=0.5)
