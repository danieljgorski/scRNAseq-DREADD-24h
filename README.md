# scRNAseq-DREADD-24h


This repository contains single-cell RNA sequencing analysis of a project led by [Dr. Katharina Botterman](mailto:katharina.bottermann@hhu.de). Samples were prepared from cardiac tissue, 24 h after ischemia/reperfusion injury from adipose specific (iAdipoQCre) hM4Di-DREADD mice. Data were processed with the Seurat toolkit, using SCTransform normalization and reference-based integration with reciprocal PCA. Control samples were used as reference, DREADD samples as query.

## Sequencing data
Sequencing data, including fastq files and count matrices will be available upon publication or request.

## Analysis
To recreate the full analysis you can follow the steps below. If you would like to process the data with your own custom workflow, a final list of cells (barcodes + metadata) after preprocessing, doublet removal and low-quality cluster removal can be found in:

* `data/basic_annotation.csv` (all interstitial cells)

### Libraries
The following scRNA-seq-specific libraries were used:

* [Seurat](https://satijalab.org/seurat/index.html) v4.1.1
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) v2.0.3
* [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) v2.13.1
* [miloR](https://marionilab.github.io/miloR/) v.1.2.0
* [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) v1.16.0
* [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) v1.22.0
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) v3.36.0
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) v4.2.1
* [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) v1.14.1
* [org.Mm.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) v3.14.0
* [yulab.utils](https://cran.r-project.org/package=yulab.utils) v0.0.4

They can be installed with the following commands:
```R
# Seurat
remotes::install_version("Seurat", version = "4.1.1")

# DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# ComplexHeatmap
library(devtools)
install_github("jokergoo/ComplexHeatmap")

# miloR
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("miloR")

# SingleCellExperiment
BiocManager::install("SingleCellExperiment")

# scater
BiocManager::install("scater")

# edgeR
BiocManager::install("edgeR")

# clusterProfiler etc.
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("org.Mm.eg.db")
install.packages("yulab.utils")
```

Additionally, the following standard R libraries were used:
* dplyr
* ggplot2
* patchwork
* readr
* ggrepel
* stringr
* knitr
* kableExtra
* forcats

### Instructions
To reproduce the analysis, clone this repository and place the cellranger aggr output inside the `data` folder. It should then contain the following:
```
scRNAseq-DREADD-24h/data/
    555_cellranger_aggr
    basic_annotation.csv
    canonical_markers.csv
    gene_signatures.xlsx
    sample-metadata.xlsx
    seurat_cell_cycle.csv
```
By starting your R session with the R project file, `scRNAseq-DREADD-24h.Rproj`, your working directory will be set to the main `scRNAseq-DREADD-24h` folder, no matter the location on your machine. This will allow easy reading/writing of data/results using relative paths.

`0-full-analysis.R` will create all necessary directories and run the full analysis in the appropriate order. Each analysis step can also be run individually for better interactivity, starting from `1-preprocessing.R`. The estimated total memory needed to hold and process the resulting Seurat objects is ~80GB.

## Examples
<p align="center">
  <img src="/examples/DimPlot_basic_annotation.png" width="1000">
</p>

<p align="center">
  <img src="/examples/Heatmap.png" width="1000">
</p>


## To-do
* Add abstract at submission
* Add full author list at submission

