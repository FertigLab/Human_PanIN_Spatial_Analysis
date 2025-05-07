---
title: "IMCpipeline_humanpanin_withfunctions_tls_analysis"
output: html_document
---

  rm(list = ls())

####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile=NULL,
                      panelDataFile=NULL,
                      dataDirectory=NULL,
                      shape_conditions=NULL,
                      color_conditions=NULL){
  ## This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ## Directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ## Read-in metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$file_name <- factor(md$file_name)
  md$File_order <- factor(md$File_order)
  md$Patient <- factor(md$Patient)
  md$Acquisition_number <- factor(md$Acquisition_number)
  md$Condition <- factor(md$Condition)
  md$Region <- factor(md$Region)
  md$TLS_status <- factor(md$TLS_status)
  md$KRAS_mut <- factor(md$KRAS_mut)
  md$Type <- factor(md$Type)
  
  rownames(md) = md$sample_id
  md$sample_id <- md$sample_id
  
  ## Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  
  ## Read fcs into fcs_raw
  fcs_raw <- read.flowSet(paste0(dataDirectory,"/",md$file_name), transformation = FALSE, truncate_max_range = FALSE)
  panel <- read_xlsx(panelDataFile)
  head(data.frame(panel))
  panel$Parameter <- gsub('-', '_', panel$Parameter)
  
  
  ## Export out the parameter/panel data from the flowFrame to edit
  ## use panel$Antigen to fix description in panel_fcs
  ## use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc)))  
  
  rownames(panel_fcs) = panel_fcs$name
  
  ## Replace desc with revised Name
  panel_fcs[panel$Parameter,]$desc<-panel$Name
  
  ## Replace parameter data in flowSet with edits
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  
  ## Assign objects to marker lists
  subtype_markers <- panel$Name[panel$Subtype == 1]
  functional_markers <- panel$Name[panel$Functional == 1]
  otherparameters <- panel$Name[panel$Other ==1]
  cluster_by <- panel$Name[panel$Cluster == 1]
  
  ## Check marker lists
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc')}
  if(!all(otherparameters %in% panel_fcs$desc)){stop('ERR: Not all otherparameters in panel_fcs$desc')}
  if(!all(cluster_by %in% panel_fcs$desc)){stop('ERR: Not all cluster markers in panel_fcs$desc')}
  
  fcs <- fsApply(fcs_raw, function(x, cofactor = 0.8){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    exprRaw<-exprs(fcs_raw[[i]])
    
    colnames(exprRaw)<-panel_fcs$desc
    
    expr<-cbind(exprs(fcs[[i]])[, union(subtype_markers,functional_markers)],exprRaw[,otherparameters])
    
    ## Combine other (spatial) data with the protein data
    colnames(expr)<-c(colnames(exprs(fcs[[i]])),otherparameters)
    
    ## Filter out any event that is 95th percentile for BOTH CD29 and CD45 (antibody aggregates)
    ##expr<-expr[expr[,"CD29"] < quantile(expr[,"CD29"], probs=0.95) & expr[,"CD45"] < quantile(expr[,"CD45"], probs=0.95),]
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs1<-flowSet(sapply(exprTr_list,flowFrame))
  
  ## Change parameter rownames
  panel_fcs1 <- pData(parameters(fcs1[[1]]))
  rownames(pData(parameters(fcs1[[1]]))) <- rownames(panel_fcs[panel_fcs$desc %in% pData(parameters(fcs1[[1]]))$name,])
  
  
  
  ###to scale every flowframe
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    expr<-exprs(fcs[[i]])
    
    expr<-t(scale(t(expr)))
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs2<-flowSet(sapply(exprTr_list,flowFrame))
  
  
  
  
  ## Get sample ids
  sample_ids <- rep(md$sample_id, fsApply(fcs1, nrow))
  
  ## Return: 
  ## fcs (only marker expressions arcsin transformed), 
  ## fcs1 (arcsin transformed + spatial parameters), 
  ## fcs2 (scaled arcsin expr per flowframe)
  ## fcsraw (all raw data), and all marker/parameter lists
  return(list('fcs'=fcs,'fcs1'=fcs1,'fcs2'=fcs2,'fcsraw'=fcs_raw,'subtype_markers'=subtype_markers,'functional_markers'=functional_markers,'otherparameters'=otherparameters,'cluster_by'=cluster_by,'sample_ids'=sample_ids,'meta_data'=md))
}


# From Melissa - I changed numclusters to 50 since that is what is currently on my HumanPanIN_merge_tls
clusterfcs <- function(fcs=output$fcs,
                       cluster_by = output$cluster_by,
                       seed=1234,plottitle='consensus_plots',
                       scaleoption=F,
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = scaleoption) %>% BuildSOM(colsToUse = cluster_by)
  
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####DIAGNOSTIC HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters,
                                            number_clusters, cluster_merging = NULL, 
                                            cluster_by=output$cluster_by,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 
  color_clusters = kovesi.rainbow_bgyrm_35_85_c69(number_clusters)
  
  
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,cluster_by]
  
  
  ## Calculate the mean expression##################################################
  
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cluster_merging$new_cluster))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cluster_merging$new_cluster,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=cluster_merging$new_cluster,
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop, 
                      gp = gpar(fill=color_list_byoriginal, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat, name="scaled",
               col=rev(brewer.rdbu(100)),
               row_order = cluster_merging[order(cluster_merging$new_cluster),]$original_cluster,
               cluster_columns = T,
               cluster_rows = T,
               border = NA,
               rect_gp = gpar(col = "white", lwd = .5),
               right_annotation = cp,
               show_row_names = T,
               row_names_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=10),
               heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
               width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=8, height=6) 
  
  return(q)
  
  dev.off() 
  
}


plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters='auto',
                                             colorbar=rev(brewer.rdbu(100)),
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = colorbar, 
                cluster_cols = T,
                cluster_rows = T, 
                labels_row = labels_row,
                #scale="row",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "white",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}



####UMAP####
#separate UMAP also created in Giotto
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~/Config/HumanPanIN_merge.xlsx',
                    seed = 1234, ncells=2000,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
  #dtf$B<- dtf$type!="B"#add a column that indicates whether the cell is a B cell or not; TRUE is non-B
  ##Why exclude B CELLS?
  ## WE HAVE NO B CELLS bc dtf$B depends on dtf$type depends on cellclustering1m which is just numbers 1:40 so..?
  #should it be the named parts cluster in merge file corresponding to it like if 30 -> grepl(cluster_merging[30,3],pattern='^B')??
  #Im blocking this out till we know why we have to do this
  ## Create subset columns without B cells (ONLY to generate the correct # of columns in inds2 object)
  #sample_ids2 <- dtf$ids[dtf$type!="B"] #sampleids without B cells
  #cell_clustering1m2 <- dtf$type[dtf$type!="B"] #clusters without B cells
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  #inds2 <- split(1:length(sample_ids2), sample_ids2) #to fill in the original indexes that do not have B cells
  samplenames <- names(inds) #create a name vector of the files
  #FOR THIS TO WORK MUST BE IN FORMAT PRE/POST Tx
  # for (i in 1:(length(samplenames)/2)){#1:15 was here because ids in dtf was 30 factors and could be divided into 0 and 6wk for each so changed it
  #   templength <- length(inds2[[i]])
  #   inds2[[i]] <- inds[[i]][dtf$B[dtf$ids==samplenames[i]]] #subsets the "B cell or not a B cell" column for each sample by the name
  #   inds2[[i]] <- inds2[[i]][1:templength]
  # }
  
  custom.settings = umap.defaults
  custom.settings$seed = seed
  
  #custom.settings$n.neighbors = neighbors
  ####umapindex generation####
  #umap ncells = table of sample ids with how many to downsample to by default col = id, row = ncells
  #sample ids = chr [1:2851129] "4927_0wk" "4927_0wk" "4927_0wk" "4927_0wk" ...
  #^ from ## Generate sample IDs corresponding to each cell in the 'expr' matrix sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  #can subset sample_ids and rerun umap 
  #can do timepoint or patient number or pre/post if you know corresp sample ids
  #sample_subset = '02_0wk' or c('02_0wk','2973_0wk') for example picking all 0 wk ones or use regex sample_ids[(grepl(sample_ids,pattern = '0wk'))]
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  
  #exclude any unassigned cluster post umap if needed--this has to be done by looking at the two columns 
  #metadata$sampleid is just a number in this metadatafile so to make unique ones combine samp_id col with timepoint (0wk)
  #to get in format of umapres2d$sample_id which looks like "02_0wk" do:
  return(umapRes2D)
}


plotUmap <- function(umapRes,seed=1234,neighbors=10,midpoint,color_clusters=colorassigned,code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    #geom_point(size = 1) +
    geom_text(aes(label = cell_clustering), size=2)+
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  
  print(ggp + facet_wrap(~ Condition, ncol=8)+ggtitle('Condition'))
  print(ggp + facet_wrap(~ Patient, ncol=5)+ggtitle('Patient'))
  print(ggp + facet_wrap(~ Region, ncol=9)+ggtitle('Region'))
  print(ggp + facet_wrap(~ TLS_status, ncol=3)+ggtitle('TLS_status'))
  
  
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()
    ) +
    
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp2)
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}





####REQUIRED LIBRARIES####
library(reshape2)
library(pals)
library(ggplot2)
library(Hmisc)
library(ComplexHeatmap)
library(ggiraphExtra)
library(diffcyt)
library(ggvoronoi)
library(ggpubr)
library(multcomp)
library(sf)
library(clusterSim)
library(circlize)
library(RColorBrewer)
library(stringr)
library(igraph)
library(readxl)
library(dplyr)
library(packcircles)
library(gridExtra)
library(limma)
library(qgraph)
library(flowCore)
library(pheatmap)
library(matrixStats)
library(ggbreak)
library(readxl)

#======================
####RUNNING DATA####
#======================



####DATA LOADING####

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()



## Read (skip if previously ran)

output <- returnfcs(metaDataFile = paste0(work,"/Config/HumanPanIn_metadata_tls.xlsx"),
                    panelDataFile = paste0(work,"/Config/HumanPanIN_panel_tls.xlsx"),
                    dataDirectory = paste0(work,"/Data"))


output<-readRDS('backup_output_40_fcs2_tls.rds')

# Make sure that your output file is using the most up-to-date metadata file you have
md <- read_excel("~/Config/HumanPanIn_metadata_tls.xlsx")
output$meta_data <- md

## Set up levels

samplevels=c("1E_25386_2_tls_ta_1",
             "1E_25386_2_tls_ta_2",
             "1M_28500_1_tls_panin_1",
             "1M_28500_1_tls_panin_2",
             "1M_28500_7_tls_panin_1",
             "2D_20447_3_tls_panin_1",
             "2HH_20447_2_tls_cpb_1",
             "2HH_20447_3_tls_cp_1",
             "2HH_20447_3_tls_cp_2",
             "2HH_20447_3_tls_cp_3",
             "2K_25590_4_tls_tb_1",
             "2O_20447_3_tls_tb_1",
             "2O_20447_4_tls_t_1",
             "2O_20447_4_tls_t_2",
             "2O_20447_5_tls_t_1",
             "2O_20447_6_tls_ta_1",
             "2O_20447_7_tls_tb_1",
             "2O_20447_8_tls_ta_1",
             "2O_20447_9_tls_ta_1",
             "2O_20447_10_tls_tb_1",
             "2Q_25590_2_tls_ta_1",
             "2Q_25590_2_tls_ta_2",
             "2Q_25590_2_tls_ta_3")

conditionlevels=c("Norm","CP", "CPB", "PanIN", "PA", "T", "TB", "TA")

patientlevels=c("1","2","3", "4", "5")

regionlevels=c("Norm","Norm_adj", "TLS", "T", "CP", "PanIN_adj", "PanIN", "T_adj", "T_bord")

tlslevels=c("NoTLS","TLSadj", "TLS")

KRASlevels=c("G12D", "G12V", "Q61H", "NA")

typelevels=c("N", "PanIN", "PDAC", "CP")


####DIAGNOSTICS####
## Spot check - number of cells per sample
cell_table <- table(output$sample_ids)
ggdf <- data.frame(sample_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$Condition <- factor(output$meta_data$Condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels = conditionlevels)
ggdf$Patient <- factor(output$meta_data$Patient[match(ggdf$sample_id,output$meta_data$sample_id)], levels= patientlevels)
ggdf$Region <- factor(output$meta_data$Region[match(ggdf$sample_id,output$meta_data$sample_id)], levels=regionlevels)
ggdf$TLS_status <- factor(output$meta_data$TLS_status[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tlslevels)
ggdf$KRAS_mut <- factor(output$meta_data$KRAS_mut[match(ggdf$sample_id,output$meta_data$sample_id)], levels=KRASlevels)
ggdf$Type <- factor(output$meta_data$Type[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)


ggp<-ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = cell_counts), angle = 45, hjust = 0.5, vjust = -0.5, size = 2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('Diagnostics_cellcounts_tls.pdf',width=20, height=20);ggp; dev.off()

## Multi-dimensional scaling plot to show similarities between samples
## Get the mean marker expression per sample
expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, fsApply(output$fcs,exprs)) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))
expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$Condition <- factor(output$meta_data$Condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels = conditionlevels)
ggdf$Patient <- factor(output$meta_data$Patient[match(ggdf$sample_id,output$meta_data$sample_id)], levels=patientlevels)
ggdf$Region <- factor(output$meta_data$Region[match(ggdf$sample_id,output$meta_data$sample_id)], levels=regionlevels)
ggdf$TLS_status <- factor(output$meta_data$TLS_status[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tlslevels)
ggdf$KRAS_mut <- factor(output$meta_data$KRAS_mut[match(ggdf$sample_id,output$meta_data$sample_id)], levels=KRASlevels)
ggdf$Type <- factor(output$meta_data$Type[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)


ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = Type, shape = Patient)) +
  geom_point(size = 2.5) +
  #geom_text(aes(label = patient_id)) +
  theme_classic()
pdf('Diagnostics_MDS_sample_tls.pdf',width=6, height=6);ggp; dev.off()

## Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
pdf('Diagnostics_Heatmap_tls.pdf',width=8, height=8)
pheatmap(expr_mean_sample_tbl[,c(2:36,39)], color = color, display_numbers = TRUE,
         number_color = "black", fontsize_number = 3, 
         clustering_method = "average")
dev.off()


####CLUSTERING####

##Revised loading depending on the diagnostics if needed

# output <- returnfcs(metaDataFile = paste0(work,"/Config/J19113_metadata.xlsx"),
# panelDataFile = paste0(work,"/Config/J19113_panel.xlsx"),
# dataDirectory = paste0(work,"/Data"))


##Clustering

output[(length(output)+1):(length(output)+3)] <- clusterfcs(fcs=output$fcs2, numclusters=40, scaleoption = T) 
#output$fcs uses just arcsin transformed data; scaleoption scales the entire dataset
#output$fcs2 uses scaled expression data that is scaled for each flowFrame

names(output)[(length(output)-2):(length(output))] <- c('code_clustering','cell_clustering','metaclusters')


####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####

## Load merge file
## Assign colors
clusterMergeFile = paste0(work,"/Config/HumanPanIN_merge_tls.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c(1:40)

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))

clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

## Metacluster heatmaps
plot_clustering_heatmap_wrapper(fcs=output$fcs2,
                                number_clusters = length(unique(output$cell_clustering1m)),
                                cell_clustering = output$cell_clustering, 
                                cluster_by = output$cluster_by,
                                clusterMergeFile = clusterMergeFile,
                                fileName = 'Clusteringheatmap_40_fcs2_tls.pdf'); dev.off()


# For merge for 40 clusters

clusterMergeFile = paste0(work,"/Config/HumanPanIN_merge_tls_fcs2_40.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("CD3- CD16+",
                "DCSIGN+",
                "CD4+ CD45RO+",
                "CD68+",
                "CD4+ FOXP3+ PD1-",
                "CD8+ CD45RO+",
                "CD8+",
                "PDPN+",
                "CD8+ GZMB+",
                "CD4+ CD45RO+ CD57+",
                "CD20+ CD45RA+",
                "CD4+ KI67+",
                "CD20+ CD21+ CD23+",
                "CD4+",
                "CD4+ CD45RO+ TOX2+ PD1+",
                "CD4+ CD45RO+ TOX2+ PD1+ KI67+",
                "CD20+ CD45RA+ KI67+",
                "Unassigned")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

plot_clustering_heatmap_wrapper2(fcs=output$fcs2,
                                 colorbar = kovesi.diverging_bwr_40_95_c42(100),
                                 subtype_markers = output$cluster_by,
                                 color_clusters = colorassigned,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 fileName = 'Clusteringheatmap_merge_40_fcs2_merge_tls.pdf');dev.off()


## Compute umap
umapRes <- do_umap(fcs=output$fcs2,
                   subtype_markers = output$subtype_markers,
                   sample_ids = output$sample_ids,
                   cell_clustering = output$cell_clustering, 
                   metadata=output$metadata,
                   clusterMergeFile=clusterMergeFile,
                   seed = 1234, 
                   ncells=500,
                   sample_subset=NULL)

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data"]]$sample_id))
umapRes$sample_id <- factor(output[["meta_data"]]$sample_id[mm], levels=samplevels)
umapRes$Patient <- factor(output[["meta_data"]]$Patient[mm], levels = patientlevels)
umapRes$Condition <- factor(output[["meta_data"]]$Condition[mm], levels=conditionlevels)
umapRes$Region <- factor(output[["meta_data"]]$Region[mm], levels=regionlevels)
umapRes$TLS_status <- factor(output[["meta_data"]]$TLS_status[mm], levels=tlslevels)
umapRes$KRAS_mut <- factor(output[["meta_data"]]$KRAS_mut[mm], levels=KRASlevels)
umapRes$Type <- factor(output[["meta_data"]]$Type[mm], levels=typelevels)



dev.off()
pdf('Umaps_50_fcs2_tls.pdf',width=20,height=20)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = colorassigned[names(colorassigned)!="NA"],
         subtype_markers = output$subtype_markers)
dev.off()



## Save output list
# saveRDS(output, file="backup_output_40_fcs2_tls.rds")
# saveRDS(umapRes, file="backup_umap_50_fcs2_tls.rds")
# umapRes<-readRDS('backup_umap_50_fcs2_tls.rds')


#plot box plots

# Subsetting TLS samples
# Looking at all cells and cell types only in TLS samples (also not intra-tumoral samples)


##Extract to get sampleIDs as a dataframe
subsetting_metadata <- as.data.frame(output$sample_ids)
##Extract TLS names from sampleID that don't have _t_
tme_cells_tls <- grep("tls", output$sample_ids,value=TRUE)
tme_cells_tls <- grep("_t_", tme_cells_tls,value=TRUE,invert=T)
# Remove this sample: 2K_25590_1_tls_1 (not actually a TLS associated with any particular region and I didn't want "NA" to show up in the plot)
tme_cells_tls <- grep("2K_25590_1_tls_1", tme_cells_tls,value=TRUE,invert=T)
# Remove this sample, TLS is too far from the PanIN
tme_cells_tls <- grep("1O_25590_6_tls_pa_1", tme_cells_tls,value=TRUE,invert=T)
# Need to remove these samples from the analysis since it turns out this patient has cholangio and not PDAC
tme_cells_tls <- grep("1D_28500_1_tls_ta_1", tme_cells_tls,value=TRUE,invert=T)
tme_cells_tls <- grep("1D_28500_1_tls_ta_2", tme_cells_tls,value=TRUE,invert=T)
tme_cells_tls <- grep("1D_28500_1_tls_ta_3", tme_cells_tls,value=TRUE,invert=T)
tme_cells_tls <- grep("1D_28500_1_tls_ta_4", tme_cells_tls,value=TRUE,invert=T)


##If the sampleID is present in the TLS names, then we write TRUE, otherwise FALSE
subsetting_metadata$include <- ifelse(subsetting_metadata$`output$sample_ids` %in% tme_cells_tls, TRUE, FALSE)
##Subset everything based on this ^ true and false vector
tme_clusters <- output$cell_clustering1m[subsetting_metadata$include]
tme_sampleids <- output$sample_ids[subsetting_metadata$include]


counts_table <- table(tme_clusters, tme_sampleids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100

counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(counts_table,file="~/counts_table_tls_fcs2_40.csv")


#Create vector of all possible celltypes
celltypes <- c("CD20+ CD21+ CD23+ KI67+", "CD20+ CD21+ CD23+", "CD20+ CD45RA+", "CD4+ CD45RO+ CCR7+ HLADR+",
               "CD8+ GZMB+", "CD68+ CD16+ CD11c+", "CD4+ CD45RO+", "CD8+ CD45RO+", "CD4+ CD45RO+ TOX2+ PD1+",
               "DCSIGN+", "CD4+ CD45RA+", "CD4+ CD45RA+ KI67+", "CD8+", "CD4+ FOXP3+ PD1-", "CD68+ CD16+ HLADR+",
               "CK+", "CD8+ GZMB+ KI67+ LAG3+", "PDPN+", "COL+ SMA+ VIM+", "Unassigned", "CD57+", "CD45+")

# specify comparisons to plot
type_comp <- list( c("N", "PanIN"), c("N", "PDAC"), c("N", "CP"), c("PanIN", "PDAC"), c("PanIN", "CP"), c("PDAC", "CP"))
type_comp_tls <- list(c("PanIN", "PDAC"), c("PanIN", "CP"), c("PDAC", "CP"))
# tls_status_comp <- 
region_comp <- list(c("PanIN", "PA"), c("T", "TB"), c("TB", "TA"), c("T", "TA"), c("T", "PanIN"))
KRAS_comp <- list(c("G12D", "G12V"), c("G12D", "Q61H"), c("G12V", "Q61H"), c("G12D", "NA"), c("G12V", "NA"), c("Q61H", "NA"))
cond_comp <- list(c("TA", "TB"), c("TB", "T"), c("T", "TA"))
cond_comp2 <- list(c("TA", "T"))
TLS_comp <- list(c("NoTLS", "TLSadj"))

# Variables in PlotPerArea function
#For each celltype, iterate through plotting
##Celltype = which cell you want to look at
##XAXIS = which column from ggdf table you would like on the xaxis
##COMPARISON = which comparison you would like stats on

PlotPerArea <- function(CELLTYPE, XAXIS, COMPARISON) {
  
  props <- as.data.frame.matrix(props[CELLTYPE,])
  counts <- as.data.frame.matrix(counts[CELLTYPE,])
  title <- CELLTYPE
  
  # Colors for this plot
  
  clusternames<-sort(unique(rownames(props)))
  
  colorassigned<-sample(kovesi.rainbow_bgyrm_35_85_c69(length(unique(rownames(props)))))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  # Create the table to plot
  ggdf <- melt(data.frame(cluster = clusternames, counts, check.names = FALSE),
               id.vars = "cluster", value.name = "counts", 
               variable.name = "sample_id")
  ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
  ggdf$Patient <- factor(output$meta_data$Patient[match(ggdf$sample_id,output$meta_data$sample_id)], levels = patientlevels)
  ggdf$Condition <- factor(output$meta_data$Condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels=conditionlevels)
  ggdf$Region <- factor(output$meta_data$Region[match(ggdf$sample_id,output$meta_data$sample_id)], levels=regionlevels)
  ggdf$TLS_status <- factor(output$meta_data$TLS_status[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tlslevels)
  ggdf$KRAS_mut <- factor(output$meta_data$KRAS_mut[match(ggdf$sample_id,output$meta_data$sample_id)], levels=KRASlevels)
  ggdf$Type <- factor(output$meta_data$Type[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
  ggdf$Area <- factor(output$meta_data$Area[match(ggdf$sample_id,output$meta_data$sample_id)])
  
  # Find the number of cells per mm2 of each region
  density <- (as.numeric(ggdf[,"counts"])*1000000)/as.numeric(as.character(ggdf[,"Area"]))
  
  density[!is.finite(density)] <- 0
  
  ggdf$density <- density
  
  # Plot in a boxplot
  
  print (ggplot(ggdf,aes(x=get(XAXIS),y=density,fill=get(XAXIS)))+
           geom_boxplot(outlier.shape=NA, lwd=0.5, color="black")+
           geom_jitter(aes(size = 12), width=0.2, size=1)+
           ggtitle(paste(title))+
           ylab("cells/mm^2")+
           # scale_y_continuous(trans = "log10")+
           # scale_y_break(c(100,200), scales = 0.5)+
           theme(axis.text.x = element_text(size=12, color="black", angle = 45, vjust = 1, hjust = 1),
                 axis.text.y = element_text(size=12, color="black"),
                 axis.title.x = element_blank(),
                 axis.line.x = element_line(size=0.5, color="black"),
                 axis.line.y = element_line(size=0.5, color="black"),
                 axis.ticks.x = element_line(size=0.5, color="black"),
                 axis.title.y = element_text(size=12, angle = 90, color="black"),
                 strip.background = element_rect(fill=NA),
                 strip.text = element_text(size=10, color="black"),
                 plot.title = element_text(size = 8, face = "bold"),
                 panel.background = element_rect(fill="white"),
                 legend.position = "none")+ 
           stat_compare_means(comparisons = COMPARISON, method = "wilcox.test", tip.length = 0,
                              label = "p.signif",
                              bracket.size = 0.5,
                              symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                 symbols = c("****", "***", "**", "*", "ns")), size = 5))
  
}


# Variables in PlotPerArea function
#For each celltype, iterate through plotting
##CELLTYPE = Celltype in the numerator
##CELLTYPE2 = Celltype in the denominator 
##XAXIS = which column from ggdf table you would like on the xaxis
##COMPARISON = which comparison you would like stats on

PlotRatioPerArea <- function(CELLTYPE, CELLTYPE2, XAXIS, COMPARISON) {
  
  props1 <- as.data.frame.matrix(props[CELLTYPE,])
  counts1 <- as.data.frame.matrix(counts[CELLTYPE,])
  props2 <- as.data.frame.matrix(props[CELLTYPE2,])
  counts2 <- as.data.frame.matrix(counts[CELLTYPE2,])
  
  title <- paste0(CELLTYPE,"/",CELLTYPE2)
  
  # Colors for this plot
  
  clusternames<-sort(unique(rownames(props1)))
  
  colorassigned<-sample(kovesi.rainbow_bgyrm_35_85_c69(length(unique(rownames(props1)))))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  # Create the table to plot for CELLTYPE
  ggdf <- melt(data.frame(cluster = clusternames, counts1, check.names = FALSE),
               id.vars = "cluster", value.name = "counts", 
               variable.name = "sample_id")
  ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
  ggdf$Patient <- factor(output$meta_data$Patient[match(ggdf$sample_id,output$meta_data$sample_id)], levels = patientlevels)
  ggdf$Condition <- factor(output$meta_data$Condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels=conditionlevels)
  ggdf$Region <- factor(output$meta_data$Region[match(ggdf$sample_id,output$meta_data$sample_id)], levels=regionlevels)
  ggdf$TLS_status <- factor(output$meta_data$TLS_status[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tlslevels)
  ggdf$KRAS_mut <- factor(output$meta_data$KRAS_mut[match(ggdf$sample_id,output$meta_data$sample_id)], levels=KRASlevels)
  ggdf$Type <- factor(output$meta_data$Type[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
  ggdf$Area <- factor(output$meta_data$Area[match(ggdf$sample_id,output$meta_data$sample_id)])
  
  # Find the number of cells per um2 of each region
  density <- (as.numeric(ggdf[,"counts"])*1000000)/as.numeric(as.character(ggdf[,"Area"]))
  
  density[!is.finite(density)] <- 0
  
  ggdf$density <- density
  
  # Create the table to plot for CELLTYPE2
  ggdf2 <- melt(data.frame(cluster = clusternames, counts2, check.names = FALSE),
                id.vars = "cluster", value.name = "counts", 
                variable.name = "sample_id")
  ggdf2$sample_id <- factor(ggdf2$sample_id, levels=samplevels)
  ggdf2$Patient <- factor(output$meta_data$Patient[match(ggdf2$sample_id,output$meta_data$sample_id)], levels = patientlevels)
  ggdf2$Condition <- factor(output$meta_data$Condition[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=conditionlevels)
  ggdf2$Region <- factor(output$meta_data$Region[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=regionlevels)
  ggdf2$TLS_status <- factor(output$meta_data$TLS_status[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=tlslevels)
  ggdf2$KRAS_mut <- factor(output$meta_data$KRAS_mut[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=KRASlevels)
  ggdf2$Type <- factor(output$meta_data$Type[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=typelevels)
  ggdf2$Area <- factor(output$meta_data$Area[match(ggdf2$sample_id,output$meta_data$sample_id)])
  
  # Find the number of cells per um2 of each region
  density <- (as.numeric(ggdf2[,"counts"])*100)/as.numeric(ggdf2[,"Area"])
  
  density[!is.finite(density)] <- 0
  
  ggdf2$density <- density
  
  # Calculate the ratio of both different cell types
  ratio <- (as.numeric(ggdf[,"density"])/as.numeric(ggdf2[,"density"]))
  
  ratio[!is.finite(ratio)] <- 0
  
  ggdf$ratio <- ratio
  
  
  
  # Plot in a boxplot
  
  print (ggplot(ggdf,aes(x=get(XAXIS),y=ratio,fill=get(XAXIS)))+
           geom_boxplot(outlier.shape=NA, lwd=0.5, color="black")+
           geom_jitter(aes(size = 12), width=0.2, size=1)+
           ggtitle(paste(title))+
           ylab("Ratio")+
           # scale_y_continuous(trans = "log10")+
           # scale_y_break(c(100,200), scales = 0.5)+
           theme(axis.text.x = element_text(size=12, color="black", angle = 45, vjust = 1, hjust = 1),
                 axis.text.y = element_text(size=12, color="black"),
                 axis.title.x = element_blank(),
                 axis.line.x = element_line(size=0.5, color="black"),
                 axis.line.y = element_line(size=0.5, color="black"),
                 axis.ticks.x = element_line(size=0.5, color="black"),
                 axis.title.y = element_text(size=12, angle = 90, color="black"),
                 strip.background = element_rect(fill=NA),
                 strip.text = element_text(size=10, color="black"),
                 plot.title = element_text(size = 8, face = "bold"),
                 panel.background = element_rect(fill="white"),
                 legend.position = "none")+ 
           stat_compare_means(comparisons = COMPARISON, method = "wilcox.test", tip.length = 0,
                              label = "p.signif",
                              bracket.size = 0.5,
                              symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                 symbols = c("****", "***", "**", "*", "ns")), size = 5))
  
}

### Cell types from merge_tls_fcs2_40 to use to make the plots running the next function
# "CD20+ CD21+ CD23+ KI67+", "CD20+ CD21+ CD23+", "CD20+ CD45RA+", "CD4+ CD45RO+ CCR7+ HLADR+",
# "CD8+ GZMB+", "CD68+ CD16+ CD11c+", "CD4+ CD45RO+", "CD8+ CD45RO+", "CD4+ CD45RO+ TOX2+ PD1+",
# "DCSIGN+", "CD4+ CD45RA+", "CD4+ CD45RA+ KI67+", "CD8+", "CD4+ FOXP3+ PD1-", "CD68+ CD16+ HLADR+",
# "CK+", "CD8+ GZMB+ KI67+ LAG3+", "PDPN+", "COL+ SMA+ VIM+", "Unassigned", "CD57+", "CD45+"

pdf('CD20+ CD45RA+_boxplots_perarea_byregion_tls_pdacptsonly.pdf',width=2.5,height=4)

PlotPerArea("CD20+ CD45RA+","Type", type_comp_tls)

dev.off()


# Look at ratio of densities

pdf('cytotoxictcells_tregs_ratio_boxplots_perarea_byregion_tls.pdf',width=2.5,height=4)

PlotRatioPerArea("CD8+ GZMB+", "CD4+ FOXP3+ PD1-","Type", type_comp_tls)

dev.off()

##### Stats table ####

# To see a table of the stats to get an exact pvalues
# Need to make the ggdf as before the function was written for making the plots


clusternames<-sort(unique(rownames(props)))

colorassigned<-sample(kovesi.rainbow_bgyrm_35_85_c69(length(unique(rownames(props)))))

names(colorassigned)<-clusternames

color_list = list(clusters=colorassigned)

# Create the table to plot
ggdf <- melt(data.frame(cluster = clusternames, counts, check.names = FALSE),
             id.vars = "cluster", value.name = "counts", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$Patient <- factor(output$meta_data$Patient[match(ggdf$sample_id,output$meta_data$sample_id)], levels = patientlevels)
ggdf$Condition <- factor(output$meta_data$Condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels=conditionlevels)
ggdf$Region <- factor(output$meta_data$Region[match(ggdf$sample_id,output$meta_data$sample_id)], levels=regionlevels)
ggdf$TLS_status <- factor(output$meta_data$TLS_status[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tlslevels)
ggdf$KRAS_mut <- factor(output$meta_data$KRAS_mut[match(ggdf$sample_id,output$meta_data$sample_id)], levels=KRASlevels)
ggdf$Type <- factor(output$meta_data$Type[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$Area <- factor(output$meta_data$Area[match(ggdf$sample_id,output$meta_data$sample_id)])

# Find the number of cells per mm2 of each region
density <- (as.numeric(ggdf[,"counts"])*1000000)/as.numeric(as.character(ggdf[,"Area"]))

density[!is.finite(density)] <- 0

ggdf$density <- density

wc_type_stats <- compare_means(density ~ Type, ggdf, method = "wilcox.test", group.by = "cluster")

