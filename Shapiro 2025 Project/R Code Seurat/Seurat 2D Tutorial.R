# General Notes 
# 
# Link to tutorial and data for visium data: https://satijalab.org/seurat/articles/spatial_vignette
# Link to CosMX Human Lung Tutorial: https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager
# 
# Human Lung dataset was analyzed with a 960 genes (assummed mRNA). There are a total of 8 samples present. Even though 960 genes were looked for on average a cell only expresses ~250. 
# This creates uncertainty. Good thing to note that in CosMX data not every cell will be expressing every gene. 

# # Loads in dataset. Only give filepath not the file
# nano.obj <- LoadNanostring(data.dir = "C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\Seurat Tutorial\\Lung5_Rep1+SMI+Flat+data\\Lung5_Rep1\\Lung5_Rep1-Flat_files_and_images", fov = "lung5_rep1")
# nano.obj <- LoadNanostring(data.dir = "/brahms/hartmana/vignette_data/nanostring/lung5_rep1", fov = "lung5.rep1")
# 
# # add in precomputed Azimuth annotations. 
# azimuth.data <- readRDS("C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\Seurat Tutorial\\nanostring_data.Rds")
# nano.obj <- AddMetaData(nano.obj, metadata = azimuth.data$annotations)
# nano.obj[["proj.umap"]] <- azimuth.data$umap
# Idents(nano.obj) <- nano.obj$predicted.annotation.l1
# 
# # set to avoid error exceeding max allowed size of globals
# options(future.globals.maxSize = 8000 * 1024^2)
# nano.obj <- SCTransform(nano.obj, assay = "Nanostring", clip.range = c(-10, 10), verbose = FALSE)
# 
# # text display of annotations and prediction scores
# head(slot(object = nano.obj, name = "meta.data")[2:5])
# 
# #C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\Seurat Tutorial\\Lung5_Rep1+SMI+Flat+data\\Lung5_Rep1\\Lung5_Rep1-Flat_files_and_images
# #C:\\Users\\ejohns\\Documents\\Shapiro Data Files\\Seurat Tutorial\\Lung5_Rep1+SMI+Flat+data
# 
# C:\Users\ejohns\Documents\Shapiro Data Files\Seurat Tutorial\Lung5_Rep1+SMI+Flat+data\Lung5_Rep1\Lung5_Rep1-Flat_files_and_images


library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
devtools::install_github('satijalab/seurat-data')
library(SeuratData)

InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot1

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))