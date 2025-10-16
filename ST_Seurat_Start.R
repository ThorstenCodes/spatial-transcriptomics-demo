# First install renv allowing to create virtual environments in R
install.packages("renv")
# Initate new environment 
renv::init()
# or activate a existing one
renv::activate()
# Now install your packages in the renv (don't use install.packages as it installs in the home directory and not in renv (not tracked then))
renv::install(c("Seurat", "tidyverse", "data.table", "lobstr"))
# Installing via BiocManager worked like this:
renv::install("BiocManager")
BiocManager::install("glmGamPoi")
renv::snapshot() # updates renv.lock file to inlcude the package


# Install some other packages
install.packages("Seurat")
install.packages("tidyverse")
install.packages("lobstr")
install.packages("data.table")
#install.packages("igraph")
#install.packages("rlang")
#install.packages("devtools")
#install.packages("BiocManager")
#BiocManager::install(c('SingleCellExperiment', 'Summarized Experiment', 'sparseMatrixStats'))
#install.packages("lsa")

install.packages('BiocManager')
BiocManager::install('glmGamPoi')

################################################################################
################################# Libraries ####################################
################################################################################

library(Seurat)
library(tidyverse)
library(data.table)
library(lobstr)
library(BiocManager)
library(glmGamPoi)


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

################################################################################
############################ Data Preprocessing ################################
################################################################################

# The Seurat RDS file is very big and takes always a long time to load and work with. 
# We will remove some large files in the next part, so that we can work better with it. 

# How big is transcript coordinates files
# We use lobstr for this as it does not require so much memory as the R internal object_size (this takes a while)

obj_size(object@misc$transcriptCoords)

# extract and safe transcript coordinates (to safe memory)
transcripts = object@misc$transcriptCoords
class(transcripts)

fwrite(transcripts, 'data/transcript_locations.csv.gz')
# clean environment
rm(transcripts)
gc(full=TRUE)
# remove transcript info from our object (since we have it saved now)
object@misc$transcriptCoords = list(NULL)
gc(full=TRUE)

# Safe object as Seurat file without transcripts (more managable because smaller)
saveRDS(object, 'data/HFC_no_transcripts.rds')

# Remove RNA_normalization from object
Assays(object)
object@assays$RNA_normalized = NULL

# We also remove the content (UMAP and PCA) from reductions as we will generate it again with our normalization
Reductions(object)
object@reductions = list(NULL)

gc(full=TRUE)
# safe again our processed seurat object
saveRDS(object, 'data/HFC_reduced.rds')

################################################################################
############################ Normalization #####################################
################################################################################
 
# If necessary reload reduced object file
object = readRDS("data/HFC_reduced.rds")


# First remove falseprobes and negativeprb from RNA file
counts = GetAssayData(object, assay = 'RNA')

# check negativeprobes
dim(counts)
# --> rows 6278 genes and 188686 cells
head(row.names(counts)) # gene names
# there rows are named negPrb (lookup in object@RNA$negprobes$data$Dimnames)
row.names(counts) %>% grep('NegPrb', ., value=TRUE)
# no ngeative probes are in our experiment, yeeeahhhh !

# check FalseCodes in the same way
row.names(counts) %>% grep('FalseCode', ., value=TRUE)
# ---> no Falsecodes inside, very good!

# No we will transform with SCTransform 

object = SCTransform(object, assay = 'RNA', new.assay.name = 'SCT')
