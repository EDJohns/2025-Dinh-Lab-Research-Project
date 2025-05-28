install.packages("remotes")
install.packages("Rtools")
install.packages("SeuratDesk")
library(Seurat)

load("C:/Users/evanj/Downloads/CosMX_protein_annot_cluster.rda")
SeuratObject <- get("CosMX_protein_annot_cluster")

loaded_names <- load("C:/Users/evanj/Downloads/CosMX_protein_annot_cluster.rda")
print(loaded_names)

library(SeuratDisk)

SaveH5Seurat(SeuratObject, filename = "SeuratObject.h5Seurat")
Convert("SeuratObject.h5Seurat", dest = "h5ad", overwrite = TRUE)

# print values
print("End of Script")


