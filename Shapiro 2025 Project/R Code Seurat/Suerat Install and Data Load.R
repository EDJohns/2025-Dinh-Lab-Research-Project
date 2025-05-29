library(Seurat)
library(SeuratObject)
library(SeuratDisk)

# Load the Seurat object from the .rda file
load("C:/Users/evanj/Downloads/CosMX_protein_annot_cluster.rda")

# Only run this if you're sure the object is named "CosMX_protein_annot_cluster"
SeuratObject <- get("CosMx_Protein_log_norm")

# Checks to see if the name loaded matches
# loaded_names <- load("C:/Users/evanj/Downloads/CosMX_protein_annot_cluster.rda")
# print(loaded_names)

# Save as .h5Seurat
SaveH5Seurat(SeuratObject, filename = "SeuratObject.h5Seurat")

# Convert to .h5ad (AnnData format for Python use)
Convert("SeuratObject.h5Seurat", dest = "h5ad", overwrite = TRUE)


SeuratObject <- DietSeurat(SeuratObject,
                           counts = TRUE,
                           data = TRUE,
                           scale.data = FALSE,
                           features = NULL,
                           assays = Assays(SeuratObject),
                           dimreducs = c("pca", "umap"),
                           graphs = NULL)


seurat_to_adata <- function(out, seu){
  seu <- Seurat::UpdateSeuratObject(seu)
  SeuratDisk::SaveH5Seurat(seu, filename = out)
  SeuratDisk::Convert(out, dest = 'h5ad')
}

seurat_to_adata("output.h5seurat", SeuratObject)