# 14_A_Summarize_Module_Scores

# summarize module score across all segments and patients

library(Seurat)
library(ggplot2)
library(ggpubr)
library(colorRamps)

sessionInfo()

plot_save_pdf <- function(plot, filename){
  ggsave(plot = plot,
         filename = paste0(filename),
         width = unit(6, "in"), height = unit(6, "in"))
}

plot_module_score_comparison <- function(seurat, subject_name, module_score, group.by,
                                         comparisons,
                                         show = TRUE, figure_path = NULL){
  cell_meta_data <- seurat@meta.data
  matlab_pal <- colorRamps::matlab.like2(length(unique(cell_meta_data[[group.by]])))
  viol <- ggplot(data = cell_meta_data, aes(x = .data[[group.by]],
                                            y = .data[[module_score]],
                                            fill = .data[[group.by]])) +
    geom_jitter(size = 0.5, color = "#000000", width = 0.05) +
    geom_violin(alpha = 0.6) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test") +
    scale_fill_manual(values = matlab_pal) +
    ggtitle(label = paste0(subject_name, " : ", module_score))
  if(show){ print(viol) }
  if(!is.null(figure_path)){
    plot_save_pdf(plot = viol,
                  filename = paste0(figure_path, "/", subject_name, "_",
                                    module_score, "_viol_", group.by, ".pdf"))
  }
}

# save module score summary info across group.by for a seurat object
module_score_summary <- function(seurat, module_score, group.by){
  mean_score <- c()
  median_score <- c()
  sd_score <- c()
  min_score <- c()
  max_score <- c()
  
  groups <- unique(seurat@meta.data[[group.by]])
  for(g in groups){
    # module score in cells belonging to group g
    group_score <- seurat@meta.data[seurat@meta.data[[group.by]] == g,][[module_score]]
    
    avg <- mean(group_score)
    med <- median(group_score)
    std <- sd(group_score)
    min_s <- min(group_score)
    max_s <- max(group_score)
    
    mean_score <- c(mean_score, avg)
    median_score <- c(median_score, med)
    sd_score <- c(sd_score, std)
    min_score <- c(min_score, min_s)
    max_score <- c(max_score, max_s)
  }
  
  df <- data.frame(
    group = groups,
    module_score = rep(module_score, length(groups)),
    mean_score = mean_score,
    median_score = median_score,
    sd_score = sd_score,
    min_score = min_score,
    max_score = max_score
  )
  
  # rename group column to the group.by argument
  colnames(df)[1] <- group.by
  
  return(df)
}

## Create Results Directories
result_dir <- "processed_data/12_A_Summarize_Module_Scores"
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

figure_dir <- "figures/12_A_Summarize_Module_Scores"
if(!dir.exists(figure_dir)){
  dir.create(figure_dir)
}

## load data
seurat <- readRDS("processed_data/09_Immune_Module_Scores/all_segments_module_scores.rds")

# create an identifier for patient number
seurat$patient <- gsub("*._", "subject_", seurat$segment)

## calcuate and save summary statitics

chemo_segment <- module_score_summary(seurat,
                                      module_score = "chemokine_module_score1",
                                      group.by = "segment")
chemo_patient <- module_score_summary(seurat,
                                      module_score = "chemokine_module_score1",
                                      group.by = "patient")

## write csvs with summary results
write.csv(chemo_segment,
          file = paste0(result_dir, "/module_score_chemokine_segment_summary.csv"))
write.csv(chemo_patient,
          file = paste0(result_dir, "/module_score_chemokine_patient_summary.csv"))

