# Volcano plots of differential gene expression analysis

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
source("scripts/etc/VolcanoPlot.R")

# Load dge data
dge_no_threshold <-
  read.csv(file = "results/12-differential-gene-expression/dge_no_threshold.csv")

# Set up output dirs
output_dirs <- c("results",
                 "results/13-volcano-plots")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Loop through clusters and output VolcanoPlots
for (i in unique(dge_no_threshold$cluster)) {
  pdf(
    file = paste0(
      "results/13-volcano-plots/VolcanoPlot_",
      i,
      ".pdf"
    ),
    height = 6,
    width = 8,
    useDingbats = F
  )
  VolcanoPlot(
    df = dge_no_threshold,
    identity = i,
    title = paste0("Differential gene expression: ", i, " - DREADD/Control")
  )
  dev.off()
}
