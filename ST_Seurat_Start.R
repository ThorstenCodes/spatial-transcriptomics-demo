# First install renv allowing to create virtual environments in R
install.packages("renv")
# Initate the environment 
renv::init()

# Install some other packages
install.packages("Seurat")
install.packages("tidyverse")
#install.packages("igraph")
#install.packages("rlang")
#install.packages("devtools")
#install.packages("BiocManager")
#BiocManager::install(c('SingleCellExperiment', 'Summarized Experiment', 'sparseMatrixStats'))
#install.packages("lsa")

################################################################################
################################# Libraries ####################################
################################################################################

library(Seurat)
library(tidyverse)


################################################################################
################################# Load Data ####################################
################################################################################

object = readRDS("data/SeuratObj_withTranscripts (1).RDS")

################################################################################
################################# Exploration ##################################
################################################################################

object@meta.data

