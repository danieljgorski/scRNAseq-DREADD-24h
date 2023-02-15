# Move project packages to HPC for manual installation

## 1. Download all project specific tarballs

  This R code will download all tarballs for the project specific pacakages
  and their dependencies, for manual installation on the HPC. It assumes the project in question has renv implemented already, which produces a project specific library found under `<project>/renv/library'.

  ```r
  library(tools)
  library(utils)
  library(BiocManager)
  project_lib <- .libPaths()[1]
  project_packages <- installed.packages(lib.loc = project_lib)
  deps <- tools::package_dependencies(packages = project_packages,
                                      which = c("Depends", "Imports"),
                                      recursive = TRUE)
  all_packs <- union(unlist(deps), project_packages)
  repos = c(BiocManager::repositories(),getOption("repos"))
  order <- download.packages(pkgs = all_packs,
                              destdir = "renv/cellar",
                              type = "source",
                              repos = repos)
  order <- basename(order[,2])
  save(order, file = "renv/cellar/order.Rdata")
  # Write a package index file in case you cannot use renv::restore on the HPC
  write_PACKAGES("renv/cellar/")
  ```

## 2. Manually move 'renv/cellar' to HPC with globus

## 3. Set the cellar path in your .Renviron file on the HPC

ssh in, move to the snakemake-node, start Miniconda/3, and activate the r4.2 conda environment, then run:

```bash
(r4.2) [gorskid@hilbert211 ~]$ echo 'RENV_PATHS_CELLAR="/gpfs/project/gorskid/scRNAseq-DREADD-24h/renv/cellar"' >> $HOME/.Renviron
```

## 4. Start an R session in the project folder

Your project folder should already be tracked by renv, you might have to to let the initial renv boostrap timeout:

```bash
(r4.2) [gorskid@hilbert211 ~]$ cd /gpfs/project/gorskid/scRNAseq-DREADD-24h/
(r4.2) [gorskid@hilbert211 scRNAseq-DREADD-24h]$ R
```

Then load `renv`, and check the cellar path:

```r
> library(renv)
> renv:::renv_paths_cellar()
[1] "/gpfs/project/gorskid/scRNAseq-DREADD-24h/renv/cellar"
```

Now restore should use the cellar tarballs first before installing from remote repos (which are not available on the HPC)

```r
> renv::restore()
```

## 5. If renv::restore() fails, proceed with install.packages()

```r
> .libPaths()
> setwd("/gpfs/project/gorskid/scRNAseq-DREADD-24h/renv/cellar")
> load("order.Rdata")
> order
> install.packages(order, repo=NULL)
```

If it turns out the dependencies for certain packages fail to install directly in R, try to `conda install r-<dependency>`. This worked for many packages e.g. `r-ggrastr`, `r-systemfonts`, `r-cairo` etc...
