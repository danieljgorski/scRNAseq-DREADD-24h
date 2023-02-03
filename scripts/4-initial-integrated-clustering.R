# Reference based integration with RPCA and clustering after doublet removal

# Load libraries
library(Seurat)
library(dplyr)

# Load object
load("results/0-objects/obj_db_removed.Rdata")

# Split object on samples
obj_list <- SplitObject(obj, split.by = "sample")

# SCTransform objects individually with sped up glmGamPoi method
obj_list <- lapply(X = obj_list, FUN = SCTransform, method = "glmGamPoi")

# Find integration features, remove DREADD feature (transgenic expression)
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3001)
features <- features[features != "R26-LSL-Gi-DREADD-genomic-sequence"]

# Prep for SCT integration, run PCA on individual samples
obj_list <- PrepSCTIntegration(object.list = obj_list,
                               anchor.features = features)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features)

# Find integration anchors, using control samples as reference
controls <- c(4, 5, 7, 8) # indices of control samples
obj_anchors <- FindIntegrationAnchors(object.list = obj_list,
                                      normalization.method = "SCT",
                                      reference = controls,
                                      anchor.features = features,
                                      dims = 1:18,
                                      reduction = "rpca",
                                      k.anchor = 5)
remove(obj)
remove(obj_list)

# Integrate data
obj <- IntegrateData(anchorset = obj_anchors,
                     normalization.method = "SCT",
                     dims = 1:18)

# Run PCA, UMAP and cluster integrated object
obj <- RunPCA(obj, npcs = 50, verbose = T)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:18, verbose = T)
obj <- FindNeighbors(obj, dims = 1:18, verbose = T)
obj <- FindClusters(obj, resolution = 0.8, verbose = T)

# Normalizing and scaling RNA assay
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     verbose = T)
obj <- ScaleData(obj, features = rownames(obj), verbose = T)

# Factor genotype level
obj@meta.data$genotype <- factor(obj@meta.data$genotype,
                                 levels = c("Control",
                                            "DREADD"))

# Save object
save(obj, file = "results/0-objects/obj_integrated.Rdata")

# Clear memory
rm(list = ls())
gc()
