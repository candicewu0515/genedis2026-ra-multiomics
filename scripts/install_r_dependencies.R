# Tested with R 4.5.2; exact package inventory is in env/R_packages.tsv.
if (getRversion() < "4.5.0") stop("Use R >=4.5; tested version 4.5.2")
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_version("Seurat", version="5.5.0", repos="https://cloud.r-project.org", upgrade="never")
if (!requireNamespace("Matrix", quietly=TRUE)) install.packages("Matrix", repos="https://cloud.r-project.org")
