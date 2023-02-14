#' Get package dependencies
#'
#' @param packs A string vector of package names
#'
#' @return A string vector with packs plus the names of any dependencies
getDependencies <- function(packs){
  dependencyNames <- unlist(
    tools::package_dependencies(packages = packs, db = installed.packages(), 
                                which = c("Depends", "Imports"),
                                recursive = TRUE))
  packageNames <- union(packs, dependencyNames)
  packageNames
}
packages <- getDependencies(rownames(installed.packages()))
repos = c(BiocManager::repositories(),getOption("repos"))
download.packages(pkgs = packages,
                             destdir = "renv/cellar",
                             type = "source",
                             repos = repos)


library(tools)
write_PACKAGES("renv/cellar/")
