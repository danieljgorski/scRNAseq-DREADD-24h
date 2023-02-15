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
write_PACKAGES("renv/cellar/")
