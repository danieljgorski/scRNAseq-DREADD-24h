# Over-representation analysis of differentially expressed genes

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(yulab.utils)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(knitr)
library(kableExtra)
source("scripts/etc/ORA_up.R")
source("scripts/etc/ORA_down.R")

# Load dge data
dge_no_threshold <-
  read.csv(file = "results/12-differential-gene-expression/dge_no_threshold.csv")

# Set up output dirs
output_dirs <- c("results",
                 "results/14-ora")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load object for background genes
load("results/0-objects/obj_annotated.Rdata")
detected <- rownames(obj)
remove(obj)

# Loop through clusters and export upregulated ORA analyses
for (i in unique(dge_no_threshold$cluster)) {
  pdf(file = paste0("results/14-ora/ORA_up_", i, ".pdf"),
      height = 4,
      width = 7,
      useDingbats = F)
  ORA_up(dge_no_threshold, i, detected, "Upregulated in DREADD")
  dev.off()
}

# Loop through clusters and export downregulated ORA analyses
for (i in unique(dge_no_threshold$cluster)) {
  pdf(file = paste0("results/14-ora/ORA_down_", i, ".pdf"),
      height = 4,
      width = 7,
      useDingbats = F)
  ORA_down(dge_no_threshold, i, detected, "Downregulated in DREADD")
  dev.off()
}
