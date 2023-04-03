# scRNAseq-DREADD-24h

This repository contains single-cell RNA sequencing analysis of a project led by [Dr. Katharina Botterman](mailto:katharina.bottermann@hhu.de). Samples were prepared from cardiac tissue, 24 h after ischemia/reperfusion injury from adipose specific (iAdipoQCre) hM4Di-DREADD mice. Data were processed with the Seurat toolkit, using SCTransform normalization and reference-based integration with reciprocal PCA. Control samples were used as reference, DREADD samples as query.

## Sequencing data

Processed and raw data have been deposited at ArrayExpress under accession E-MTAB-12882, a public link will be available upon publication or request.

## Analysis

To recreate the full analysis you can follow the steps below. If you would like to process the data with your own custom workflow, a final list of cells (barcodes + metadata) after preprocessing, doublet removal and low-quality cluster removal can be found in:

* `data/basic_annotation.csv`

### R & Libraries

R version 4.1.2 was used, packages versions and external sources were recorded with `renv` v0.16.0. To create a library with matching package versions from the `renv.lock` file, start an R session in the project directory and run:

```r
renv::restore()
```

More information on ```renv``` can be found [here](https://rstudio.github.io/renv/articles/renv.html).

### Data

To reproduce the analysis, place the extracted `cellranger aggr` output inside the `data` folder. It should then contain the following:

```bash
scRNAseq-DREADD-24h/data/
    555_cellranger_aggr
    basic_annotation.csv
    canonical_markers.csv
    gene_signatures.xlsx
    sample-metadata.xlsx
    seurat_cell_cycle.csv
```

By starting your R session with the R project file, `scRNAseq-DREADD-24h.Rproj`, your working directory will be set to the main `scRNAseq-DREADD-24h` folder, no matter the location on your machine. This will allow easy reading/writing of data/results using relative paths.

`0-full-analysis.R` will create all necessary directories and run the full analysis in the appropriate order. Each analysis step can also be run individually for better interactivity, starting from `1-preprocessing.R`. The estimated total memory needed to hold and process the resulting Seurat objects is ~140GB.

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
