# Basic cluster annotation based on canonical cell-type markers

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/8-basic-annotation")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

# Load object
load("results/0-objects/obj_integrated_clean.Rdata")

# Calculate seurat_cluster markers, with high stringency to return very
# specific markers, quickly.
seurat_cluster_markers <- FindAllMarkers(obj,
                                  assay = "RNA",
                                  logfc.threshold = 1.5,
                                  min.pct = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.0001,
                                  densify = T,
                                  verbose = T)
write.csv(seurat_cluster_markers,
          file = "results/8-basic-annotation/seurat_cluster_markers.csv",
          row.names = F)
top5 <- seurat_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Read in canonical markers
markers <- read.csv(file = "data/canonical_markers.csv")
markers <- markers$All

# Loop through markers and generate Feature and VlnPlots of expression
count = 0
total = length(markers)
for (i in markers) {
  count = count + 1
  p1 <- FeaturePlot(obj, features = i, label = T, raster = F)
  p2 <- VlnPlot(obj, features = i, pt.size = 0, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/8-basic-annotation/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1, 2)))
  dev.off()
  print(paste(i, " - done ", count, "/", total))
}

# Export a basic DimPlot of clusters
p <- DimPlot(obj, label = T, raster = F)
pdf(file = "results/8-basic-annotation/DimPlot.pdf",
    useDingbats = F)
print(p)
dev.off()

###############################################################################
# Annotation notes
###############################################################################
# 
# 7 - EC - Cdh5	Kdr	Pecam1
#
# 20 - Mural (Pericytes and SMC) - Myh11	Rgs5	Pdgfrb	Vtn	Cspg4	
#
# 10 - B-cell-1 - Cd79a Ms4a1 H2-Aa Igkc
# 22 - B-cell-2 - Cd79a Ms4a1 H2-Aa Igkc, S100a8+?, S100a9+? Csf3r+?
# 23 - B-cell-3 - Cd79a Ms4a1 H2-Aa Igkc
#
# 15 - T-cell - Cd3d, Cd3e, Lef1
#
# 6 - Fibro-Rest - Col1a1, Tcf21-high, Pdfgra, Gsn, Cthrc1-, Postn-low
# 14 - Fibro-Act - Col1a1, Tcf21-low, Pdgfra, Cthrc1-low, Postn-low, Timp1
# 19 - Fibro-Myo - Col1a1-high, Tcf21-low, Pdgfra, Cthrc1-high, Postn-high, 
#                  Mik67+ contains cycling fibroblasts aswell.
#
# 0 - Gran-1 - S100a8	S100a9	Csf3r	
# 1 - Gran-2 - S100a8	S100a9	Csf3r	
# 3 - Gran-3 - S100a8	S100a9	Csf3r	
# 9 - Gran-4 - S100a8	S100a9	Csf3r	
# 13 - Gran-5 - S100a8	S100a9	Csf3r	
#
# 12 - DC-1 - Cd209a Cd74 H2-Aa	Ccr2	Cd68	Fcgr1
# 21 - DC-2 - Cd209a Cd74 H2-Aa	Ccr2	Cd68	
#
# 8 - Mac-Mono - Mix of bone marrow derived monocytes and differentiating macs:
#               Ly6c2+ Ccr2-high Fcgr1+ Cd68+, low expression of classic 
#               mac/phagocytosis markers (Adgre1-low, H2-Aa-low, Maf-low,
#               Mertk-low, Trem2-low),low expression of complement 
#               (C1qa-low, C1qb-low, C1qc-low)
# 2 - Mac-1 - monocyte-derived pro-inflammatory classical M1 macrophages:
#               Ly6c2- Ccr2+ Fcgr1+ Cd68+ Adgre1+ H2-Aa- Cx3cr1- Maf+ Mertk+ 
#               Trem2+ Il1b+
# 4 - Mac-2 - monocyte-derived pro-inflammatory classical M1 macrophages:
#               Ly6c2- Ccr2+ Fcgr1+ Cd68+ Adgre1-low H2-Aa-high Maf+ Mertk-low 
#               Trem2-low Il1b+ these, are pro-inflammatory, but have high 
#               H2-Aa, perhaps they are the first M2 transitioners
# 5 - Mac-3 - monocyte-derived pro-inflammatory classical M1 macrophages:
#               Ly6c2- Ccr2+ Fcgr1+ Cd68+ Adgre1-low H2-Aa- Cx3cr1- Maf+ Mertk+
#               Trem2+ Il1b+
# 16 - Mac-4 - monocyte-derived pro-inflammatory macrophages:
#               Ly6c2- Ccr2+ Fcgr1+ Cd68+ Adgre1+ H2-Aa-low Cx3xr1+ Timd4-
#               Lyve1- Maf+ Mertk+ Trem2+ Il1b+ C1qb-high C1qa-high
# 18 - Mac-5 - pro-inflammatory macrophages, but not phagocytotic :
#               Ly6c2- Ccr2-low Fcgr1-low Cd68-low Adgre1-low Maf-low Mertk-
#               Trem2-low, Il1b-low, C1qa-high, C1qb-high Pf4-high
# 11 - Mac-TR - Tissue resident macrophages:
#               Ly6c2- Fcgr1+ Cd68+ Adgre1-high Ccr2- Timd4+ Lyve1+ Cx3cr1+
#               [Sub-pop is Timd4- Lyve1- H2-Aa-high like Farbehi et al. 2019])
# 17 - Mac-IFN - Interferon responsive macrophages:
#               Ly6c2+ Fcgr1+ Cd68+ Adgre1-low Ifit3+ Ifi204+ Ifi209+

###############################################################################

# Rename Idents to annotations
obj <- RenameIdents(obj,
                    "0" = "Gran-1",
                    "1" = "Gran-2",
                    "2" = "Mac-1",
                    "3" = "Gran-3",
                    "4" = "Mac-2",
                    "5" = "Mac-3",
                    "6" = "Fibro-Rest",
                    "7" = "EC",
                    "8" = "Mac-Mono",
                    "9" = "Gran-4",
                    "10" = "B-cell-1",
                    "11" = "Mac-TR",
                    "12" = "DC-1",
                    "13" = "Gran-5",
                    "14" = "Fibro-Act",
                    "15" = "T-cell",
                    "16" = "Mac-4",
                    "17" = "Mac-IFN",
                    "18" = "Mac-5",
                    "19" = "Fibro-Myo",
                    "20" = "Mural",
                    "21" = "DC-2",
                    "22" = "B-cell-2",
                    "23" = "B-cell-3")

# Store renamed idents as a new meta data column, set as Idents
obj@meta.data$basic_annotation <- Idents(obj)

# Refactor annotation levels
source("scripts/etc/dimplotlevels.R")
obj@meta.data$basic_annotation <- factor(obj@meta.data$basic_annotation,
                                         levels = dimplotlevels)
DimPlot(obj,
        group.by = "basic_annotation",
        label = T,
        repel = T)

# Set Idents as re-factored basic_annotation identities
Idents(obj) <- "basic_annotation"

# Save object with basic annotations
save(obj, file = "results/0-objects/obj_annotated.Rdata")

# Saved basic annotation, barcodes and UMAP embeddings etc. for consistency in
# external usage, de-comment if you want to overwrite 
# barcodes <- rownames(obj@meta.data)
# annotation <- obj@meta.data$basic_annotation
# genotype <- obj@meta.data$genotype
# sample <- obj@meta.data$sample
# UMAP_1 <- Embeddings(obj[["umap"]])[,1]
# UMAP_2 <- Embeddings(obj[["umap"]])[,2]
# basic_annotation <- data.frame(barcodes,
#                                annotation,
#                                genotype,
#                                sample,
#                                UMAP_1,
#                                UMAP_2)
# write.csv(basic_annotation, file = "data/basic_annotation.csv", row.names = F)
