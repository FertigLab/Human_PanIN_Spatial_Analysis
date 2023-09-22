# 16_A_TLS_Selection_for_CoGAPS

library(Seurat)
library(CoGAPS)
library(Matrix)
library(dplyr)

sessionInfo()

# create results directory
result_dir <- "processed_data/16_A_TLS_Selection_for_CoGAPS"
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

# load seurat data for each subject
ser01 <- readRDS("processed_data/15_A_Picking_Neighbor_Spots/subject_01_adjacent_spots.rds")
ser02 <- readRDS("processed_data/15_A_Picking_Neighbor_Spots/subject_02_adjacent_spots.rds")
ser03 <- readRDS("processed_data/15_A_Picking_Neighbor_Spots/subject_03_adjacent_spots.rds")
ser04 <- readRDS("processed_data/15_A_Picking_Neighbor_Spots/subject_04_adjacent_spots.rds")
ser05 <- readRDS("processed_data/15_A_Picking_Neighbor_Spots/subject_05_adjacent_spots.rds")

# merge seurat objects 
seurat <- merge(ser01, c(ser02, ser03, ser04, ser05))

# remove individual objects from memory and save the aggregated data
rm(ser01)
rm(ser02)
rm(ser03)
rm(ser04)
rm(ser05)

saveRDS(seurat, file = "processed_data/15_A_Picking_Neighbor_Spots/all_segments_adjacent_spots.rds")

# limit the spots to TLS and surroundings
seurat <- seurat[, seurat$TLS_neighbors_2 %in% c("TLS", "TLS_neighbor")]

# check seurat dimension
dim(seurat)
table(seurat$TLS_neighbors_2)

# save a plot of picked spots to confirm
Idents(seurat) <- seurat$TLS_neighbors_2
p <- SpatialDimPlot(seurat, ncol = 4)
ggsave(plot = p,
       file = paste0(result_dir, "/TLS+neighbors_for_CoGAPS.pdf"),
       width = unit(12, "in"), height = unit(30, "in"))

# save just the TLS and adjacent splots
saveRDS(seurat, file = "processed_data/15_A_Picking_Neighbor_Spots/TLS_adjacent_spots.rds")

## save log-transformed counts matrix
# extract counts matrix
mat <- seurat@assays$Spatial@counts

# Preprocessing on expression matrix for CoGAPS

# matrix features before any subsetting or transformations
dim(mat)
class(mat)
range(mat)
mean(mat)

# Remove cells with no signal (all expression values in column are 0)
mat <- mat[, apply(mat, 2, max) > 0]
dim(mat)
# Remove genes with no signal/constant signal (stdev of row is 0)
mat <- mat[apply(mat, 1, sd) > 0, ]
dim(mat)
# Remove genes with no signal (max of row is 0)
mat <- mat[apply(mat, 1, max) > 0, ]
dim(mat)
# Log transform based on the counts data being whole numbers 
mat <- log2(mat+1)

# Check min of transformed data
min(mat)
# Recheck the dim
dim(mat)
class(mat)
range(mat)
mean(mat)

# Save row and column names for use in CoGAPS parameter creation
geneNames<- rownames(mat)
sampleNames<- colnames(mat)

saveRDS(geneNames, file = paste0(result_dir,
                                 "/TLS_visium_geneNames.rds"))
saveRDS(sampleNames, file = paste0(result_dir,
                                   "/TLS_visium_sampleNames.rds"))

# Save the matrix as a csv
mat_dense <- as.matrix(mat)
rownames(mat_dense) <- geneNames
colnames(mat_dense) <- sampleNames

write.csv(
  mat_dense,
  file = paste0(result_dir, "/TLS_visium_logCounts.csv"))

# create params objects 

target_patterns <- c(4, 6, 8, 10, 12)

for(t in target_patterns){
  params <- CogapsParams(
    nPatterns = t,
    nIterations = 50000,
    seed = 1234,
    sparseOptimization = TRUE,
    distributed = "genome-wide",
    geneNames = geneNames,
    sampleNames = sampleNames
  )
  params <- setDistributedParams(params, nSets = 16)
  print(params)
  
  # save parameters
  saveRDS(params, file = 
            paste0(result_dir, "/CoGAPS_parameters_", t, "_patterns.rds"))
}

# clear environment
rm(list = ls())
