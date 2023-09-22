# 06_Identify_Pt_VarGenes
# Jacob Mitchell

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(ggrepel)

# create directories
data_dir <- "processed_data/05_Aggregate_by_Patient"
results_dir <- "processed_data/06_Identify_Pt_VarGenes"
figure_dir <- "figures/06_Identify_Pt_VarGenes"

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}
if(!dir.exists(figure_dir)){
  dir.create(figure_dir, recursive = TRUE)
}

# list seurat files with data 
seg_dir_list <- c(
  "processed_data/05_Aggregate_by_Patient/subject_01_aggregated_segments.rds",
  "processed_data/05_Aggregate_by_Patient/subject_02_aggregated_segments.rds",
  "processed_data/05_Aggregate_by_Patient/subject_03_aggregated_segments.rds",
  "processed_data/05_Aggregate_by_Patient/subject_04_aggregated_segments.rds",
  "processed_data/05_Aggregate_by_Patient/subject_05_aggregated_segments.rds"
)

for(seg_dir in seg_dir_list){
  seg_name <- gsub(paste0(data_dir, "/"), "", seg_dir)
  seg_name <- gsub("_aggregated_segments.rds", "", seg_name)
  
  seg_figure_dir <- paste0(figure_dir, "/", seg_name)
  if(!dir.exists(seg_figure_dir)){
    dir.create(seg_figure_dir)
  }
  
  print(paste0("Preprocessing: ", seg_name))
  
  ## read in file
  ser <- readRDS(seg_dir)
  
  writeLines(paste0(
    "nCount_Spatial Summary:\n",
    capture.output(summary(ser$nCount_Spatial))[1], "\n",
    capture.output(summary(ser$nCount_Spatial))[2], "\n",
    "nFeature_Spatial Summary:\n",
    capture.output(summary(ser$nFeature_Spatial))[1], "\n",
    capture.output(summary(ser$nFeature_Spatial))[2]
  ),
  con = paste0(seg_figure_dir, "/", seg_name, "_qc.txt"))
  
  # plot QC parameter measures
  ## nCount
  vln_count <- VlnPlot(ser, 
                       features = "nCount_Spatial",
                       group.by = "cell_type_confirmed",
                       pt.size = 0.1) +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1)) +
    NoLegend()
  
  ggsave(filename = paste0(seg_figure_dir, "/", seg_name, "_nCount_QC.pdf"),
         plot = vln_count,
         width = unit(8, "in"),
         height = unit(6, "in"))
  
  ## nFeature
  vln_feat <- VlnPlot(ser, 
                       features = "nFeature_Spatial",
                       group.by = "cell_type_confirmed",
                       pt.size = 0.1) +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1)) +
    NoLegend()
  
  ggsave(filename = paste0(seg_figure_dir, "/", seg_name, "_nFeature_QC.pdf"),
         plot = vln_feat,
         width = unit(8, "in"),
         height = unit(6, "in"))
  
  # renormalize across all segments from the subject
  # SCT transform
  ser <- SCTransform(ser, assay = "Spatial")
  
  ## plot normalized counts
  ser <- GroupCorrelation(ser, group.assay = "Spatial",
                          assay = "SCT", slot = "scale.data", do.plot = FALSE)
  sct_plot <- GroupCorrelationPlot(ser, assay = "SCT", 
                                   cor = "nCount_Spatial_cor") +
    ggtitle("SCTransform Normalization") +
    theme(plot.title = element_text(hjust = 0.5))
  sct_plot
  ggsave(filename = paste0(seg_figure_dir, "/", seg_name, "_SCT_QC.pdf"),
         plot = sct_plot,
         width = unit(6, "in"),
         height = unit(6, "in"))
  
  # Identify highly-variable features
  ser <- FindVariableFeatures(ser, selection.method = 'vst', nfeatures = 2000)
  
  ser <- FindSpatiallyVariableFeatures(ser, assay = "SCT",
                                       features = VariableFeatures(ser)[1:500], 
                                       selection.method = "markvariogram")
  top_features <- head(SpatiallyVariableFeatures(ser,
                                                 selection.method = "markvariogram"), 10)
  
  # Save intermediate file after data scaling
  saveRDS(ser,
          file = paste0(results_dir, "/", seg_name, "_VarFeatures.rds"))
  print("==========")
  print(paste0(seg_name, " preprocessing complete"))
  print(paste0(seg_name, " top spatially variable features: "))
  print(top_features)
  print("==========")
  print("")
}

print("All segment preprocessing complete")
