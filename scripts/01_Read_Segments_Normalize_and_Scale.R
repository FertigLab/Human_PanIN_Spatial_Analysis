# 01_Read_Segments_Normalize_and_Scale
# Jacob Mitchell

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(ggrepel)

# create directories
data_dir <- "data/MultiLesionCohort"
results_dir <- "processed_data/01_Read_Segments_Normalize_and_Scale"
figure_dir <- "figures/01_Read_Segments_Normalize_and_Scale"

if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}
if(!dir.exists(figure_dir)){
  dir.create(figure_dir, recursive = TRUE)
}

# list directories containing 10x output files
seg_dir_list <- list.dirs(data_dir, recursive = FALSE)
seg_dir_list <- seg_dir_list[grepl("CP", seg_dir_list) |
                               grepl("PANIN", seg_dir_list) |
                               grepl("PDAC", seg_dir_list) |
                               grepl("NRL", seg_dir_list)]

for(seg_dir in seg_dir_list){
  seg_name <- gsub(paste0(data_dir, "/"), "", seg_dir)
  
  seg_figure_dir <- paste0(figure_dir, "/", seg_name)
  if(!dir.exists(seg_figure_dir)){
    dir.create(seg_figure_dir)
  }
  
  print(paste0("Preprocessing: ", seg_name))
  
  ## read in file
  ser <- Load10X_Spatial(seg_dir)
  
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
  vln_count <- VlnPlot(ser, features = "nCount_Spatial", pt.size = 0.1) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    NoLegend()
  space_count <- SpatialFeaturePlot(ser, features = "nCount_Spatial",
                                    pt.size.factor = 2.0) +
    scale_fill_viridis_c() +
    theme(legend.position = "bottom",
          legend.text = element_text(angle = -45, vjust = 1, hjust = 0))
  wrap_plots(vln_count, space_count)
  ggsave(filename = paste0(seg_figure_dir, "/", seg_name, "_nCount_QC.pdf"),
         plot = wrap_plots(vln_count, space_count),
         width = unit(8, "in"),
         height = unit(6, "in"))
  
  ## nFeature/nCount
  vln_feat <- VlnPlot(ser, features = "nFeature_Spatial", pt.size = 0.1) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    NoLegend()
  sctr_feat <- FeatureScatter(ser,
                              feature1 = "nCount_Spatial",
                              feature2 = "nFeature_Spatial") +
    NoLegend()
  space_feat <- SpatialFeaturePlot(ser, features = "nFeature_Spatial",
                                   pt.size.factor = 2.0) +
    scale_fill_viridis_c() +
    theme(legend.position = "bottom",
          legend.text = element_text(angle = -45, vjust = 1, hjust = 0))
  wrap_plots(vln_feat, sctr_feat, space_feat, ncol = 2)
  ggsave(filename = paste0(seg_figure_dir, "/", seg_name, "_nFeature_QC.pdf"),
         plot = wrap_plots(vln_feat, sctr_feat, space_feat),
         width = unit(12, "in"),
         height = unit(6, "in"))
  
  # remove spots with no counts
  ser <- subset(ser, subset = nCount_Spatial > 0)
  
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
  top_var <- head(VariableFeatures(ser), 10)
  
  ser <- FindSpatiallyVariableFeatures(ser, assay = "SCT",
                                        features = VariableFeatures(ser)[1:500], 
                                        selection.method = "markvariogram")
  top_features <- head(SpatiallyVariableFeatures(ser,
                                                 selection.method = "markvariogram"), 10)
  top_feat_plot <- SpatialFeaturePlot(ser,
                                      features = top_features,
                                      ncol = 5, alpha = c(0.1,1))
  
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
