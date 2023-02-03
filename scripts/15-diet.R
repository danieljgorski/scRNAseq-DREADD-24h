# Create a smaller diet version of the full annotated object

# Load Library
library(Seurat)

# Load object
load("results/0-objects/obj_annotated.Rdata")

# Keep only necessary counts and data slots
obj <- DietSeurat(obj,
                  counts = F,
                  data = T,
                  scale.data = F,
                  assays = "RNA",
                  dimreducs = "umap")

# Save clean object
save(obj, file = "results/0-objects/obj_annotated_diet.Rdata")
