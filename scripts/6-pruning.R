# Pruning low-quality clusters after initial clustering

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(ggrepel)

# Load object
load("results/0-objects/obj_integrated.Rdata")


# Explore at risk clusters that could be aggregated low-quality cells or
# heterotypic doublets
###############################################################################
# At risk low quality clusters include:
# 8, 11, 24, 16, 22, 18, 10, 21, 25

# 8 This cluster is Pecam1+ and Kdr+, likely ECs, ECs are always a little more 
#stressed (Higher mt, lower nFeature than other intersitial cells), but these 
# are healthy ECs. Will leave in.

# 11 This cluster is Cd79a+ and Ms4a1+, likely B-Cells, B-cells always have 
# low nFeatures, they are not high expressors. Will leave in.

# 24 This cluster sits between ECs and granulocytes, has typical low nFeatures
# of granulocytes, is S100a8+ and S100a9+ but Kdr- and Pecam1-. But this has a
# markedly higher percent mt than other granulocytes (0,9,1,4) and a lower
# amount of marker genes. Will remove.

# 16 This cluster is Col1a1+, Sparc+, Bgn+ so it is likely fibroblasts, but 
# has markedly higher percent mt and lower nFeatures than other fibroblasts
# (7, 23, 14). These are likely stress fibros, will remove.

# 22 This cluster is Pecam1+ and Kdr+, likely ECs with slightly higher
# stress or different function than larger EC cluster (8). The marker genes
# are clearly EC based, but the nFeatures is also markedly lower than 8. 
# Will remove.

# 18 This cluster is S100a8+ and S100a9+, likely stressed granulocytes with
# high percent mt and very low cluster marker genes. Will remove.

# 10 This cluster has a smattering of Cd68 and Fcgr1 expression, high percent mt
# with low nFeature and markers genes (Gm26917, Gm42418, Lars2). A similar type
# of cluster has been published before, but I believe these are stressed
# macrophages. Will remove.

# 21 This cluster is Kdr+ and Pecam1+ but also is Ctla2a+, a marker of T-cells.
# Likely multiplets or stressed ECs, will remove.

# 25 This cluster is Actc1+ Tmp1+ Tnnt2+, likely cardiomyocytes which are 
# inherently stressed or cell debris from broken CMs that made it through FACS
# and droplet formation. Will remove.

# Possible heterotypic doublet/triplet clusters include:
# 20

# 20 is a combination of ECs (Kdr+, Pecam1+) and granulocytes (S100a8+, 
# S100a9+). Likely extravasating, will remove.
###############################################################################

# Exclude low-quality clusters determined above
obj <- subset(x = obj,
              idents = c("24", "16", "22", "18", "10", "21", "25", "20"),
              invert = TRUE)

# Remove extra meta data columns created during integrated clustering
obj@meta.data <- select(obj@meta.data,
                        -c(nCount_SCT,
                           nFeature_SCT,
                           integrated_snn_res.0.8,
                           seurat_clusters))

# Keep only necessary counts and data slots
obj <- DietSeurat(obj,
                  counts = T,
                  data = T,
                  scale.data = F,
                  features = NULL,
                  assays = "RNA",
                  dimreducs = NULL,
                  graphs = NULL)

# Save clean object
save(obj, file = "results/0-objects/obj_clean.Rdata")

# Clear memory
rm(list = ls())
gc()
