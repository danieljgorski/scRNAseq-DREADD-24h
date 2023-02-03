# Preprocessing

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/0-objects",
                 "results/1-preprocessing")

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
library(dplyr)

# Read in data
counts <- Read10X("data/555_cellranger_aggr/outs/count/filtered_feature_bc_matrix")
obj <- CreateSeuratObject(counts = counts, 
                          project = "555",
                          min.cells = 3, 
                          min.features = 200, 
                          names.field = 2, 
                          names.delim = "-")
remove(counts)

# Adding metadata
obj <- AddMetaData(obj, metadata = "60min", col.name = "ischemia")

obj <- AddMetaData(obj, metadata = "24h", col.name = "reperfusion")

obj@meta.data <- obj@meta.data %>%
  mutate(isolation = case_when(orig.ident %in% c("1", "2", "3", "4") ~
                                 "08_Dec_2021",
                               orig.ident %in% c("5", "6", "7", "8") ~
                                 "04_May_2022"))

obj@meta.data <- obj@meta.data %>%
  mutate(bmfz_sample = case_when(orig.ident == "1" ~ "555-1",
                                 orig.ident == "2" ~ "555-2",
                                 orig.ident == "3" ~ "555-3",
                                 orig.ident == "4" ~ "555-4",
                                 orig.ident == "5" ~ "555-9",
                                 orig.ident == "6" ~ "555-10",
                                 orig.ident == "7" ~ "555-11",
                                 orig.ident == "8" ~ "555-12"))

obj@meta.data <- obj@meta.data %>%
  mutate(sample = case_when(orig.ident == "1" ~ "555-1-hM4Di-het",
                            orig.ident == "2" ~ "555-2-hM4Di-wt",
                            orig.ident == "3" ~ "555-3-hM4Di-het",
                            orig.ident == "4" ~ "555-4-hM4Di-wt",
                            orig.ident == "5" ~ "555-9-hM4Di-wt",
                            orig.ident == "6" ~ "555-10-hM4Di-het",
                            orig.ident == "7" ~ "555-11-hM4Di-het",
                            orig.ident == "8" ~ "555-12-hM4Di-wt"))


obj@meta.data <- obj@meta.data %>%
  mutate(animal_number = case_when(orig.ident == "1" ~ "752646",
                                   orig.ident == "2" ~ "752647",
                                   orig.ident == "3" ~ "752651",
                                   orig.ident == "4" ~ "754524",
                                   orig.ident == "5" ~ "782953",
                                   orig.ident == "6" ~ "782956",
                                   orig.ident == "7" ~ "782964",
                                   orig.ident == "8" ~ "782965"))

obj@meta.data <- obj@meta.data %>%
  mutate(genotype_long = case_when(orig.ident %in% c("2", "4", "5", "8") ~
                                "Adipoq_CreERT2_het_hM4Di_flox_wt",
                               orig.ident %in% c("1", "3", "6", "7") ~
                                "Adipoq_CreERT2_het_hM4Di_flox_het"))

obj@meta.data <- obj@meta.data %>%
  mutate(genotype = case_when(orig.ident %in% c("2", "4", "5", "8") ~
                                "Control",
                              orig.ident %in% c("1", "3", "6", "7") ~
                                "DREADD"))

# Quality control filtering
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
cells_pre_qc <- length(colnames(obj))

# Pre-filter QC metrics by sample
p <- VlnPlot(obj,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "sample",
             pt.size = 0)
pdf(file = "results/1-preprocessing/VlnPlot_QC_metrics_pre-filter.pdf",
    width = 12,
    height = 6,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# Filter
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)
cells_post_qc <- length(colnames(obj))
percent_passed <- (cells_post_qc / cells_pre_qc) * 100
qc_filter <- as.data.frame(cells_pre_qc)
qc_filter$cells_post_qc <- cells_post_qc
qc_filter$percent_passed <- percent_passed
write.csv(qc_filter,
          file = "results/1-preprocessing/qc_filter.csv",
          row.names = F)

# Post-filter QC metrics by sample
p <- VlnPlot(obj,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "sample",
             pt.size = 0)
pdf(file = "results/1-preprocessing/VlnPlot_QC_metrics_post-filter.pdf",
    width = 12,
    height = 6,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# Save object
save(obj, file = "results/0-objects/obj_preprocessed.Rdata")

# Clear memory
rm(list = ls())
gc()
