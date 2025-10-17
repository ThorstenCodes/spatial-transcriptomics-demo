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
# The data set is too big to run memory-intensive SCT-Transformation on the Mac Book. Thus, we will only run a subset of cells and fovs!
# Furthermore we will remove some large files which we then will calculate again (for demonstration).


# The next line is only to show how big the S4 Seurat object is!
# We use lobstr for this as it does not require so much memory as the R internal object_size (this takes a while)

obj_size(object@misc$transcriptCoords)

### Tasks ###
# 1.) Determine how many slides there are and how many fovs are on each slide using the metadata
# 2.) We want to work with around 10 000 cells. How many Fovs we need to include to have enough cells. Do the fov have similar cell numbers?
# 3.) Subset the object to the reduced cell number (see next Task section for more details)

# Solution to 1.)

number_of_slides <- object@meta.data %>% distinct(slide_ID_numeric) %>% count()
number_of_slides

# --> There is only one slide! Good let's check how many FOVs are on the slide

metadata %>% filter(slide_ID_numeric == 1) %>% distinct(fov) %>%  count() 
#392

# 2.) How many of the FOVs can we use to have 10K cells?
fovs <- metadata %>% pull(fov) %>% unique()

# initialize empty list to store results
cells_per_fov <- list()  

for (f in fovs) {
  subset_metadata <- metadata %>% 
    filter(fov == f)
  cells_per_fov[[as.character(f)]] <- nrow(subset_metadata)
}

class(cells_per_fov)
fov_counts <- unlist(cells_per_fov)
# cumulative sum of cell
cumsum(fov_counts)

# Find how many FOVs you can subtract before exceeding 10000
# This will generate a vector from the point where the argument is true but we want only to see the first number in this vector, thus [1]
which(cumsum(fov_counts) > 10000)[1] - 1

# so 22 FOVs are below 10000 and 23 cross the threshold slightly.

# 3.) Subsetting the data
### Tasks ####
# 1.) Isolate cell_ids for all cells detected in the selected Fovs. 
# 2.) Filter the count matrix for the cell ids 
# 3.) SCT Transform the smaller dataset on your Computer

#Filter the name of all cells from the first 22 FOVs (100-122)
cells_subset <- metadata %>% select(cell_id, fov) %>% filter(fov %in% (100:122) ) %>% pull(cell_id)
class(cells_subset)

# Subset cells from object using Seurats internal subset function (does it on all values!)
object <- subset(object, cells = cells_subset)
View(object_sub)

# Finally I remove some data from the object as I will calculate it again and to make the object smaller.
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
object@reductions = list()

gc(full=TRUE)
# safe again our processed seurat object
saveRDS(object, 'data/HFC_reduced_subset.rds')


################################################################################
############################ Normalization #####################################
################################################################################
 
# First remove false codes and negative probes from RNA file
counts = GetAssayData(object, assay = 'RNA')

# check negative probes
dim(counts)
# --> rows 6278 genes and 188686 cells
head(row.names(counts)) # gene names
# there rows are named negPrb (look up in object@RNA$negprobes$data$Dimnames)
row.names(counts) %>% grep('NegPrb', ., value=TRUE)
# RESULT: no negative probes are in our experiment!

# check False Codes in the same way
row.names(counts) %>% grep('FalseCode', ., value=TRUE)
# ---> no False codes inside, very good!


# Now we will transform with SCTransform 
object <-  SCTransform(object, assay = 'RNA', new.assay.name = 'SCT', )

### Tasks ###
# 1.) Now run the PCA
# 2.) Illustarte how the PCA in a coordinate system
# 3.) Run UMAP

# Run PCA on subset
object <-  RunPCA(object, assay = 'SCT', reduction.name = 'PCA', npcs = 50) # npcs 50 is default, but one could run less

# Plot PCA
DimPlot(object, reduction = 'PCA')

# Run and Plot UMAP
object <- RunUMAP(object, reduction = 'PCA', reduction.name = 'UMAP', dims = 1:30, repulsion.strength = 5)

DimPlot(object, reduction = 'UMAP')


