# Quality control after initial clustering

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/5-post-clustering-qc")

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
library(ggplot2)
library(patchwork)
library(readr)
library(ggrepel)

# Load object
load("results/0-objects/obj_integrated.Rdata")

# Switch default assay to RNA
DefaultAssay(obj) <- "RNA"

# Dissociation related genes
diss_genes <- unique(c("Atf3", "Btg2", "Cebpb", "Cebpb",
                       "Cxcl3", "Cxcl2", "Cxcl1",
                       "Dnaja1", "Dnajb1", "Dusp1",
                       "Egr1", "Fos", "Fosb", "Hsp90aa1",
                       "Hsp90ab1", "Hspa1a", "Hspa1b",
                       "Hspa1a", "Hspa1b", "Hspa8",
                       "Hspb1", "Hspe1", "Hsph1", "Id3",
                       "Ier2", "Jun", "Junb", "Jund",
                       "Mt1", "Nfkbia", "Nr4a1", "Ppp1r15a",
                       "Socs3", "Zfp36"))
obj <- AddModuleScore(obj,
                      features = list(diss_genes),
                      ctrl = 50,
                      name = "diss_genes")
p1 <- FeaturePlot(obj, features = "diss_genes1", label = T, raster = F) +
  ggtitle("Dissociation gene score")
p2 <- VlnPlot(obj, features = "diss_genes1", pt.size = 0, sort = T) +
  ggtitle("Dissociation gene score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/5-post-clustering-qc/dissocation_genes.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# Batch effect exploration
p <- DimPlot(obj, group.by = "sample", raster = F)
pdf(file = "results/5-post-clustering-qc/DimPlot_sample.pdf")
print(p)
dev.off()

p <- DimPlot(obj, group.by = "isolation", raster = F)
pdf(file = "results/5-post-clustering-qc/DimPlot_isolation.pdf")
print(p)
dev.off()

p <- DimPlot(obj, group.by = "genotype", raster = F)
pdf(file = "results/5-post-clustering-qc/DimPlot_genotype.pdf")
print(p)
dev.off()

# nFeature_RNA
p1 <- FeaturePlot(obj, features = "nFeature_RNA", label = T, raster = F)
p2 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/5-post-clustering-qc/nFeature_RNA.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# nCount_RNA
p1 <- FeaturePlot(obj, features = "nCount_RNA", label = T, raster = F)
p2 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/5-post-clustering-qc/nCount_RNA.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# percent.mt
p1 <- FeaturePlot(obj, features = "percent.mt", label = T, raster = F)
p2 <- VlnPlot(obj, features = "percent.mt", pt.size = 0, sort = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NoLegend()
pdf(file = "results/5-post-clustering-qc/percent_mt.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

# DREADD expression check
p <- VlnPlot(obj,
             "R26-LSL-Gi-DREADD-genomic-sequence",
             group.by = "sample",
             sort = T)
pdf(file = "results/5-post-clustering-qc/dreadd_by_sample.pdf",
    width = 14,
    height = 6,
    useDingbats = F)
print(p)
dev.off()

# Cell-cycle phase
cc_genes <- read_csv("data/seurat_cell_cycle.csv")
s_genes <- as.character(cc_genes$mouse.s.genes)
g2m_genes <- as.character(cc_genes$mouse.g2m.genes)
remove(cc_genes)
obj <- CellCycleScoring(obj,
                        s.features = s_genes,
                        g2m.features = g2m_genes)
obj@meta.data$Phase <- factor(obj@meta.data$Phase,
                              levels = c("G1", "S", "G2M"))
phase_membership <- (prop.table(table(obj$Phase, Idents(obj)),
                                margin = 2) * 100)
phase_membership <- as.data.frame(phase_membership)
colnames(phase_membership) <- c("Phase", "Cluster", "Percent")
phase_membership$percent_round <- round(phase_membership$Percent)
write.csv(phase_membership,
          file = "results/5-post-clustering-qc/phase_membership.csv",
          row.names = F)
phase_membership$Phase <- factor(phase_membership$Phase,
                                 levels = c("G1", "S", "G2M"))
p1 <- DimPlot(obj,
              group.by = "Phase",
              cols = c("#e5e5e5", "#3a86ff", "#ffaa00"),
              raster = F)
p2 <- ggplot(phase_membership,
             aes(fill = Phase,
                 y = Percent,
                 x = Cluster)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Percent of total") +
  xlab("Identity") +
  geom_text(aes(label = paste0(percent_round, "%", sep = "")),
            position = position_stack(vjust = 0.5), size = 2.5) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "Black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = c("#e5e5e5", "#3a86ff", "#ffaa00"))
pdf(file = "results/5-post-clustering-qc/cell_cycle_phase.pdf",
    width = 17,
    height = 6,
    useDingbats = F)
p1 + p2 + plot_layout(guides = "collect")
dev.off()

# Find markers of initial clustering
cluster_markers <- FindAllMarkers(obj,
                                  assay = "RNA",
                                  logfc.threshold = 0.5,
                                  min.pct = 0.5,
                                  only.pos = T,
                                  return.thresh = 0.001,
                                  densify = T)
write.csv(cluster_markers,
          file = "results/5-post-clustering-qc/cluster_markers.csv",
          row.names = F)
top5 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,
          file = "results/5-post-clustering-qc/cluster_markers_top_5.csv",
          row.names = F)
top10 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top5,
          file = "results/5-post-clustering-qc/cluster_markers_top_10.csv",
          row.names = F)

# Counting cluster markers
marker_count <- cluster_markers %>%
  group_by(cluster) %>%
  summarise(nMarkers = n())
marker_count$cluster <- as.character(marker_count$cluster)

# Summarizing quality control
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
qc_summary <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(mean(percent.mt),
            mean(diss_genes1),
            mean(nFeature_RNA),
            mean(nCount_RNA),
            calculate_mode(Phase))
qc_summary <- qc_summary %>% left_join(marker_count,
                                      by = c("seurat_clusters" = "cluster"))
colnames(qc_summary) <- c("cluster",
                          "mean_percent_mt",
                          "mean_dissociation_score",
                          "mean_nFeatures",
                          "mean_nCounts",
                          "cell_cycle",
                          "nMarkers")
write.csv(qc_summary,
          file = "results/5-post-clustering-qc/qc_summary.csv",
          row.names = F)

# Plotting QC summary
p <- ggplot(qc_summary, aes(x = mean_percent_mt,
                       y = mean_nFeatures,
                       colour = nMarkers,
                       label = cluster)) +
  scale_colour_viridis_c() +
  geom_point(aes(shape = cell_cycle), size = 3) +
  geom_text_repel() +
  ggtitle("Initial clustering QC")
pdf(file = "results/5-post-clustering-qc/qc_summary_scatter.pdf",
    width = 8,
    height = 6,
    useDingbats = F)
print(p)
dev.off()

# Plotting QC overview summary + nFeature + Clusters
p1 <- DimPlot(obj, label = T, raster = F) + NoLegend()
p2 <- FeaturePlot(obj,
                  features = "nFeature_RNA",
                  label = T,
                  raster = F) + NoLegend()
p3 <- ggplot(qc_summary, aes(x = mean_percent_mt,
                            y = mean_nFeatures,
                            colour = nMarkers,
                            label = cluster)) +
  scale_colour_viridis_c() +
  geom_point(aes(shape = cell_cycle), size = 3) +
  geom_text_repel() +
  ggtitle("Initial clustering QC")
pdf(file = "results/5-post-clustering-qc/qc_overview.pdf",
    width = 13,
    height = 10,
    useDingbats = F)
print((p1 / p2) | p3 )
dev.off()

# Simple annotation, to identify remaining possible heterotypic multiplets

# Read in markers
markers <- read.csv(file = "data/canonical_markers.csv")
markers <- markers$All

# Loop through markers and generate Feature and VlnPlots of expression
for (i in markers) {
  p1 <- FeaturePlot(obj, features = i, label = T, raster = F)
  p2 <- VlnPlot(obj, features = i, pt.size = 0, sort = T) + NoLegend() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  pdf(file = paste0("results/5-post-clustering-qc/Feat_VlnPlot_", i, ".pdf"),
      width = 16,
      height = 5,
      useDingbats = F)
  print(p1 + p2 + plot_layout(ncol = 2, widths = c(1,2)))
  dev.off()
}

# Save object
save(obj, file = "results/0-objects/obj_integrated.Rdata")

# Clear memory
rm(list = ls())
gc()
