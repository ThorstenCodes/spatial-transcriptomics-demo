#Converting .rds file to h5ad file for python on MacOS

#Install packages required for conversion:
install.packages('Seurat')
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
install.packages("reticulate")   # needed for Python bridge
install.packages("remotes")      # needed to install from GitHub
remotes::install_github("cellgeni/sceasy")
py_install("anndata")

# Load them
library(Seurat)
library(Matrix)
library(sceasy)
library(SingleCellExperiment)
library(reticulate)

setwd('Your_path_to_.rda_file') # Caution it is rda not rds !

load('colon_reference_sce.rda')

sceasy::convertFormat(
  colon_reference_sce,
  from = "sce",
  to = "anndata",
  outFile = "colon_reference.h5ad"
)
