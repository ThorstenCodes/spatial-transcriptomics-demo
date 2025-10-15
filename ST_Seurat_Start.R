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
############################ Data Exploration ##################################
################################################################################

metadata <- object@meta.data

# quick way to find how many FOV
length(unique(metadata$fov))
# brings the same result (was the video way)
slide_fov = metadata %>% 
  select(fov, slide_ID_numeric) %>% 
  distinct()
# so in total we have 392 Field of Views (FOV)!
# How much cells are in a field of view, exemplary for the 1st FOV
# if you want the number of FOV in the slide directly add "%>% count()" or one can also use "%>% nrow()"
# but since i use the variable later i need the rest of the dataframe
slide_fov1 = metadata %>% 
  filter(fov == 1,
         slide_ID_numeric == 1)
## --> output is 424, so we have 424 cells in the first FOV (prepare with all FOVs, big dataset incoming!)

# Lets make a plot of the FOV1 only (each pixel equals 0.18 µm):
# These are all centroids of all cells!
# cell types annotated
slide_fov1 %>% 
  ggplot() +
  geom_point(aes(x = x_FOV_px, y= y_FOV_px, color= RNA_nbclust_clusters_refined), size=0.5) +
  coord_fixed()

# cluster assignments 
slide_fov1 %>% 
  ggplot() +
  geom_point(aes(x = x_FOV_px, y= y_FOV_px, color= spatialClusteringAssignments), size=0.5) +
  coord_fixed()
# when we would add the polygon information we would see the cells on slide 1 FOV 1.

#################################
#### Exploring Single Cell ######
#################################


object@misc$transcriptCoords %>% 
  filter(slideID == 1,
         fov == 1) %>%  pull(CellId) %>%  unique() %>% sort() %>% head()

### Cell 1 is not in the FOV so we look at cell 2

slide1_fov1_cell2 = object@misc$transcriptCoords %>% 
  filter(slideID == 1,
         fov == 1, CellId == 2) #--> 192 transcripts!
# Plot cell 2

slide1_fov1_cell2 %>% 
  ggplot() +
  geom_point(aes(x= x_FOV_px, y= y_FOV_px, color=target), size= 1)+
  coord_fixed() +
  theme(legend.position = "none")




