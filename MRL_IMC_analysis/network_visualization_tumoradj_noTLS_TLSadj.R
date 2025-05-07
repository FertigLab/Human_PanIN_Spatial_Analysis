#####DISTANCE RELATIONSHIP ANALYSIS#####

library(qgraph)
library(spatstat)
library(stringr)
library(pals); library(RColorBrewer)
library(pheatmap); library(ComplexHeatmap)
library(scales)
library(readxl)
library(dplyr); library(reshape2)
library(ggplot2)
library(flowCore)
library(Hmisc)

rm(list = ls())

workd<-getwd()

output<-readRDS('backup_output_50_fcs2_merge4.rds')

#KEY POINTS IN THE CODE THAT TOGGLE THE EXCLUSION OF TLS REGIONS ARE MARKED ***

#####LEVELS AND RELEVANT CELL TYPES#####

clusterlevels <- c("CD8+",
                   "CD8+ CD45RO+",
                   "CD8+ GZMB+",
                   "CD8+ GZMB+ KI67+ LAG3+",
                   "CD4+ CD45RA+",
                   "CD4+ CD45RA+ KI67+",
                   "CD4+ CD45RO+ CCR7+ HLADR+",
                   "CD4+ CD45RO+",
                   "CD4+ CD45RO+ TOX2+ PD1+",
                   "CD4+ FOXP3+ PD1-",
                   "CD20+ CD45RA+",
                   "CD20+ CD21+ CD23+",
                   "CD57+",
                   "CD57+ TOX2+",
                   "DCSIGN+",
                   "CD68+",
                   "CD68+ CD16+ HLADR+",
                   "CD68+ CD16+ CD11c+",
                   "PDPN+",
                   "COL+ SMA+ VIM+",
                   "CK+",
                   "CK+ PDL1+",
                   "KI67+",
                   "CD45+",
                   "Unassigned")

unassigned<-c("Unassigned","CD45+","KI67+")

#####COLORS AND LEGENDS#####

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
colorassigned<-hue_pal()(length(clusterlevels)-length(unassigned))
names(colorassigned)<-clusterlevels[clusterlevels %nin% unassigned]
hex <- hue_pal()(9)
colorassignedbroad <- c(rep(hex[1],4), #CD8
                        rep(hex[2],5), #CD4
                        rep(hex[3],2), #FOXP3
                        rep(hex[4],5), #CD20
                        rep(hex[5],2), #CD57
                        rep(hex[6],3), #Myeloid
                        rep(hex[7],2), #CK
                        rep(hex[8],2)) #Stroma
                        
allcelltypes<-clusterlevels[clusterlevels %nin% unassigned]
legendctype<-as.data.frame(cbind(paste0("ctype",1:length(allcelltypes)),allcelltypes))
legendctype$maintype<-1
legendctype$maintype[str_detect(legendctype$allcelltypes,"CD4+")]<-"CD4+"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CD8+")]<-"CD8+"
legendctype$maintype[str_detect(legendctype$allcelltypes,"FOXP3+")]<-"FOXP3+"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CD20+")]<-"CD20+"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CD57+")]<-"CD57+"
legendctype$maintype[str_detect(legendctype$allcelltypes,"DCSIGN+")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CD68+")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"PDPN+")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"COL+")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CK+")]<-"CK"

includectypes<-str_sub(legendctype$V1,6,)

colorlist <- cbbPalette[as.numeric(as.factor(legendctype$maintype))]


#####COUNTS AND PERCENTAGES#####

counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdft <- melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
              id.vars = "cluster", value.name = "counts", 
              variable.name = "sample_id")
ggdft$sample_id <- factor(ggdft$sample_id)
ggdft$cluster <- factor(ggdft$cluster)
ggdft$Type <- factor(output$meta_data$Type[match(ggdft$sample_id,output$meta_data$sample_id)])
ggdft$TLS_status <- factor(output$meta_data$TLS_status[match(ggdft$sample_id,output$meta_data$sample_id)])
ggdft$Region <- factor(output$meta_data$Region[match(ggdft$sample_id,output$meta_data$sample_id)])


totalcounts_celltype_type <- ggdft %>% group_by(cluster, Type) %>% summarize_at(vars(counts),funs(sum))

# write.csv(totalcounts_celltype_type, paste0(workd,"/Results/Totalcounts_type.csv"))

totalcounts_celltype_tls <- ggdft %>% group_by(cluster, TLS_status) %>% summarize_at(vars(counts),funs(sum))

# write.csv(totalcounts_celltype_tls, paste0(workd,"/Results/Totalcounts_TLS.csv"))

totalcounts_celltype_region <- ggdft %>% group_by(cluster, Region) %>% summarize_at(vars(counts),funs(sum))

# write.csv(totalcounts_celltype_region, paste0(workd,"/Results/Totalcounts_region.csv"))

totalcounts<-ggdft %>% group_by(cluster, Type, TLS_status, Region) %>% summarize_at(vars(counts),funs(sum))

# write.csv(totalcounts, paste0(workd,"/Results/Totalcounts.csv"))


#***TO EXCLUDE TLS REGIONS
ggdft <- ggdft[ggdft$TLS_status!="TLS",]

# Remove all panin adjacent regions, just want to evaluate real panin regions
ggdft <- ggdft[ggdft$Region!="PanIN_adj",]

# Exclude intra tumoral regions since you want to compare to only tumor border or tumor adjacent
ggdft <- ggdft[ggdft$Region!="T",]

totalcounts<-ggdft %>% group_by(cluster, Type, TLS_status) %>% summarize_at(vars(counts),funs(sum))
totalcounts_celltype_type <- ggdft %>% group_by(cluster, Type) %>% summarize_at(vars(counts),funs(sum))
totalcounts_celltype_tls <- ggdft %>% group_by(cluster, TLS_status) %>% summarize_at(vars(counts),funs(sum))

#percentage of each respective total

totalcounts_Normal <- totalcounts_celltype_type[totalcounts_celltype_type$Type=="N",]
totalcounts_CP <- totalcounts_celltype_type[totalcounts_celltype_type$Type=="CP",]
totalcounts_PanIN <- totalcounts_celltype_type[totalcounts_celltype_type$Type=="PanIN",]
totalcounts_PDAC <- totalcounts_celltype_type[totalcounts_celltype_type$Type=="PDAC",]

totalNormal<-sum(totalcounts_Normal$counts)
totalCP<-sum(totalcounts_CP$counts)
totalPanIN<-sum(totalcounts_PanIN$counts)
totalPDAC<-sum(totalcounts_PDAC$counts)

pct_Normal <- totalcounts_Normal$counts/totalNormal*100
names(pct_Normal)<- totalcounts_Normal$cluster
pct_CP <- totalcounts_CP$counts/totalCP*100
names(pct_CP)<- totalcounts_CP$cluster
pct_PanIN <- totalcounts_PanIN$counts/totalPanIN*100
names(pct_PanIN)<- totalcounts_PanIN$cluster
pct_PDAC <- totalcounts_PDAC$counts/totalPDAC*100
names(pct_PDAC)<- totalcounts_PDAC$cluster

ggdf_pctNormal <- melt(pct_Normal);ggdf_pctNormal$Type<-"N"
ggdf_pctCP <- melt(pct_CP);ggdf_pctCP$Type<-"CP"
ggdf_pctPanIN <- melt(pct_PanIN);ggdf_pctPanIN$Type<-"PanIN"
ggdf_pctPDAC <- melt(pct_PDAC);ggdf_pctPDAC$Type<-"PDAC"

ggdf_pct<-rbind(ggdf_pctNormal,ggdf_pctCP,ggdf_pctPanIN,ggdf_pctPDAC)
ggdf_pct$cluster<-c(rownames(ggdf_pctNormal),rownames(ggdf_pctCP),rownames(ggdf_pctPanIN),rownames(ggdf_pctPDAC))
ggdf_pct$cluster<-factor(ggdf_pct$cluster, levels=clusterlevels)
rownames(ggdf_pct)<-1:nrow(ggdf_pct)
ggdf_pct$Type<-factor(ggdf_pct$Type, levels=c("N","PanIN","PDAC","CP"))

ggdf_pct_sub <- ggdf_pct[ggdf_pct$cluster %nin% unassigned,]

bp <- ggplot(ggdf_pct_sub, aes(x = Type, y = value, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  scale_fill_manual(values = hue_pal()(length(unique(ggdf_pct$cluster))))+
                    #breaks = clusterlevels,
                    #labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

pdf(paste0(workd,"/Results/Abundance_stackedbar_Type_tumoradj_noTLS.pdf"), width=7, height=5); bp; dev.off()

#identify which cell types in either of the data subsets are very rare to exclude - have less than 0.01%

include_Normal<-which(pct_Normal>0.01)
include_Normal<-legendctype$V1[match(names(include_Normal), legendctype$allcelltypes)]
include_Normal<-include_Normal[!is.na(include_Normal)]

include_CP<-which(pct_CP>0.01)
include_CP<-legendctype$V1[match(names(include_CP), legendctype$allcelltypes)]
include_CP<-include_CP[!is.na(include_CP)]

include_PanIN<-which(pct_PanIN>0.01)
include_PanIN<-legendctype$V1[match(names(include_PanIN), legendctype$allcelltypes)]
include_PanIN<-include_PanIN[!is.na(include_PanIN)]

include_PDAC<-which(pct_PDAC>0.01)
include_PDAC<-legendctype$V1[match(names(include_PDAC), legendctype$allcelltypes)]
include_PDAC<-include_PDAC[!is.na(include_PDAC)]

#get out expression levels, X, and Y coords

expr <- fsApply(output$fcs1, exprs) #create expression matrix

expr0<-data.frame(expr[,c(union(output$subtype_markers,output$functional_markers),"CellId","X_coord","Y_coord")],
                  cluster=output$cell_clustering1m,
                  sample_id=output$sample_ids)
expr0$sample_id <- factor(expr0$sample_id)
expr0$cluster <- factor(expr0$cluster)
expr0$Type <- factor(output$meta_data$Type[match(expr0$sample_id,output$meta_data$sample_id)])
expr0$TLS_status <- factor(output$meta_data$TLS_status[match(expr0$sample_id,output$meta_data$sample_id)])
expr0$Region <- factor(output$meta_data$Region[match(expr0$sample_id,output$meta_data$sample_id)])

expr0<-expr0[expr0$cluster %nin% unassigned,]

#***TO EXCLUDE TLS REGIONS
expr0<-expr0[expr0$TLS_status!= "TLS",]

# Remove all panin adjacent regions, just want to evaluate real panin regions
expr0 <- expr0[expr0$Region!="PanIN_adj",]

# Exclude intra tumoral regions since you want to compare to only tumor border or tumor adjacent
expr0 <- expr0[expr0$Region!="T",]


expr0_Normal<-expr0[expr0$Type=="N",]
expr0_CP<-expr0[expr0$Type=="CP",]
expr0_PanIN<-expr0[expr0$Type=="PanIN",]
expr0_PDAC<-expr0[expr0$Type=="PDAC",]

#####CREATE DISTANCE MATRICES FOR KEY TYPES#####

Normal<-unique(expr0[expr0$Type=="N",]$sample_id)
CP<-unique(expr0[expr0$Type=="CP",]$sample_id)
PanIN<-unique(expr0[expr0$Type=="PanIN",]$sample_id)
PDAC<-unique(expr0[expr0$Type=="PDAC",]$sample_id)

##Normal

expr_Normal<-c()

for(k in 1:length(Normal)){
  expr_k<-expr0[expr0$sample_id==Normal[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy <- matrix(nrow=nrow(expr_k),ncol=length(clusterlevels)) 
  colnames(dummy) <- paste0("ctype",1:length(clusterlevels))
  dummy <- as.data.frame(dummy)

  expr_k <- data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  distByCelltype <- distByCelltype[,clusterlevels]
  colnames(distByCelltype) <- paste0("ctype",1:length(clusterlevels))
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  
  expr_Normal <- rbind(expr_Normal, expr_k)
}

expr_Normal_m<-as.matrix(expr_Normal[,colnames(expr_Normal)[str_detect(colnames(expr_Normal),"ctype")]])
expr_Normal$ctype_no<-paste0("ctype",match(expr_Normal$cluster,legendctype$allcelltypes))
rownames(expr_Normal_m)<-expr_Normal$ctype_no
expr_Normal_m[is.infinite(expr_Normal_m)]<-NA

expr_Normal<-expr_Normal_m
expr_Normal[expr_Normal == 0] <- NA

saveRDS(expr_Normal,'backup_dist_Normal_tumoradj_noTLS.rds')

#load previously saved distance matrix

expr_Normal<- readRDS('backup_dist_Normal_tumoradj_noTLS.rds')

#create expression + distance data frames

expr_Normalcombined <- cbind(expr0_Normal,expr_Normal)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_Normalrm0<-expr_Normal[,colnames(expr_Normal) %nin% colnames(expr_Normal)[colSums(expr_Normal, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types

expr_Normalrm<-expr_Normalrm0[rownames(expr_Normalrm0) %in% include_Normal, colnames(expr_Normalrm0) %in% include_Normal]

#aggregate the distances and clean up the matrix
mat_Normal=aggregate(x=expr_Normalrm, by=list(rownames(expr_Normalrm)), FUN=mean, na.rm=T)
groupnames<-mat_Normal$Group.1
mat_Normal<-as.matrix(mat_Normal[,2:ncol((mat_Normal))])
rownames(mat_Normal)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_Normal)<-str_sub(colnames(mat_Normal),6,) #simplify colnames
mat_Normalex<-mat_Normal[rownames(mat_Normal)[c(rownames(mat_Normal) %in% includectypes)],colnames(mat_Normal)[c(colnames(mat_Normal) %in% includectypes)]] #excluding Other subtypes
dist_Normal<-mat_Normalex


##CP

expr_CP<-c()

for(k in 1:length(CP)){
  expr_k<-expr0[expr0$sample_id==CP[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy <- matrix(nrow=nrow(expr_k),ncol=length(clusterlevels)) 
  colnames(dummy) <- paste0("ctype",1:length(clusterlevels))
  dummy <- as.data.frame(dummy)
  
  expr_k <- data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  distByCelltype <- distByCelltype[,clusterlevels]
  colnames(distByCelltype) <- paste0("ctype",1:length(clusterlevels))
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  
  expr_CP <- rbind(expr_CP, expr_k)
}

expr_CP_m<-as.matrix(expr_CP[,colnames(expr_CP)[str_detect(colnames(expr_CP),"ctype")]])
expr_CP$ctype_no<-paste0("ctype",match(expr_CP$cluster,legendctype$allcelltypes))
rownames(expr_CP_m)<-expr_CP$ctype_no
expr_CP_m[is.infinite(expr_CP_m)]<-NA

expr_CP<-expr_CP_m
expr_CP[expr_CP == 0] <- NA

saveRDS(expr_CP,'backup_dist_CP_tumoradj_noTLS.rds')

#load previously saved distance matrix

expr_CP<- readRDS('backup_dist_CP_tumoradj_noTLS.rds')

#create expression + distance data frames

expr_CPcombined <- cbind(expr0_CP,expr_CP)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_CPrm0<-expr_CP[,colnames(expr_CP) %nin% colnames(expr_CP)[colSums(expr_CP, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types

expr_CPrm<-expr_CPrm0[rownames(expr_CPrm0) %in% include_CP, colnames(expr_CPrm0) %in% include_CP]

#aggregate the distances and clean up the matrix
mat_CP=aggregate(x=expr_CPrm, by=list(rownames(expr_CPrm)), FUN=mean, na.rm=T)
groupnames<-mat_CP$Group.1
mat_CP<-as.matrix(mat_CP[,2:ncol((mat_CP))])
rownames(mat_CP)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_CP)<-str_sub(colnames(mat_CP),6,) #simplify colnames
mat_CPex<-mat_CP[rownames(mat_CP)[c(rownames(mat_CP) %in% includectypes)],colnames(mat_CP)[c(colnames(mat_CP) %in% includectypes)]] #excluding Other subtypes
dist_CP<-mat_CPex


##PanIN

expr_PanIN<-c()

for(k in 1:length(PanIN)){
  expr_k<-expr0[expr0$sample_id==PanIN[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy <- matrix(nrow=nrow(expr_k),ncol=length(clusterlevels)) 
  colnames(dummy) <- paste0("ctype",1:length(clusterlevels))
  dummy <- as.data.frame(dummy)
  
  expr_k <- data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  distByCelltype <- distByCelltype[,clusterlevels]
  colnames(distByCelltype) <- paste0("ctype",1:length(clusterlevels))
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  
  expr_PanIN <- rbind(expr_PanIN, expr_k)
}

expr_PanIN_m<-as.matrix(expr_PanIN[,colnames(expr_PanIN)[str_detect(colnames(expr_PanIN),"ctype")]])
expr_PanIN$ctype_no<-paste0("ctype",match(expr_PanIN$cluster,legendctype$allcelltypes))
rownames(expr_PanIN_m)<-expr_PanIN$ctype_no
expr_PanIN_m[is.infinite(expr_PanIN_m)]<-NA

expr_PanIN<-expr_PanIN_m
expr_PanIN[expr_PanIN == 0] <- NA

saveRDS(expr_PanIN,'backup_dist_PanIN_tumoradj_noTLS.rds')

#load previously saved distance matrix

expr_PanIN<- readRDS('backup_dist_PanIN_tumoradj_noTLS.rds')

#create expression + distance data frames

expr_PanINcombined <- cbind(expr0_PanIN,expr_PanIN)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_PanINrm0<-expr_PanIN[,colnames(expr_PanIN) %nin% colnames(expr_PanIN)[colSums(expr_PanIN, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types

expr_PanINrm<-expr_PanINrm0[rownames(expr_PanINrm0) %in% include_PanIN, colnames(expr_PanINrm0) %in% include_PanIN]

#aggregate the distances and clean up the matrix
mat_PanIN=aggregate(x=expr_PanINrm, by=list(rownames(expr_PanINrm)), FUN=mean, na.rm=T)
groupnames<-mat_PanIN$Group.1
mat_PanIN<-as.matrix(mat_PanIN[,2:ncol((mat_PanIN))])
rownames(mat_PanIN)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_PanIN)<-str_sub(colnames(mat_PanIN),6,) #simplify colnames
mat_PanINex<-mat_PanIN[rownames(mat_PanIN)[c(rownames(mat_PanIN) %in% includectypes)],colnames(mat_PanIN)[c(colnames(mat_PanIN) %in% includectypes)]] #excluding Other subtypes
dist_PanIN<-mat_PanINex



##PDAC

expr_PDAC<-c()

for(k in 1:length(PDAC)){
  expr_k<-expr0[expr0$sample_id==PDAC[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy <- matrix(nrow=nrow(expr_k),ncol=length(clusterlevels)) 
  colnames(dummy) <- paste0("ctype",1:length(clusterlevels))
  dummy <- as.data.frame(dummy)
  
  expr_k <- data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  distByCelltype <- distByCelltype[,clusterlevels]
  colnames(distByCelltype) <- paste0("ctype",1:length(clusterlevels))
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  
  expr_PDAC <- rbind(expr_PDAC, expr_k)
}

expr_PDAC_m<-as.matrix(expr_PDAC[,colnames(expr_PDAC)[str_detect(colnames(expr_PDAC),"ctype")]])
expr_PDAC$ctype_no<-paste0("ctype",match(expr_PDAC$cluster,legendctype$allcelltypes))
rownames(expr_PDAC_m)<-expr_PDAC$ctype_no
expr_PDAC_m[is.infinite(expr_PDAC_m)]<-NA

expr_PDAC<-expr_PDAC_m
expr_PDAC[expr_PDAC == 0] <- NA

saveRDS(expr_PDAC,'backup_dist_PDAC_tumoradj_noTLS.rds')

#load previously saved distance matrix

expr_PDAC<- readRDS('backup_dist_PDAC_tumoradj_noTLS.rds')

#create expression + distance data frames

expr_PDACcombined <- cbind(expr0_PDAC,expr_PDAC)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_PDACrm0<-expr_PDAC[,colnames(expr_PDAC) %nin% colnames(expr_PDAC)[colSums(expr_PDAC, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types

expr_PDACrm<-expr_PDACrm0[rownames(expr_PDACrm0) %in% include_PDAC, colnames(expr_PDACrm0) %in% include_PDAC]

#aggregate the distances and clean up the matrix
mat_PDAC=aggregate(x=expr_PDACrm, by=list(rownames(expr_PDACrm)), FUN=mean, na.rm=T)
groupnames<-mat_PDAC$Group.1
mat_PDAC<-as.matrix(mat_PDAC[,2:ncol((mat_PDAC))])
rownames(mat_PDAC)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_PDAC)<-str_sub(colnames(mat_PDAC),6,) #simplify colnames
mat_PDACex<-mat_PDAC[rownames(mat_PDAC)[c(rownames(mat_PDAC) %in% includectypes)],colnames(mat_PDAC)[c(colnames(mat_PDAC) %in% includectypes)]] #excluding Other subtypes
dist_PDAC<-mat_PDACex




#####VISUALIZE NETWORKS#####

pdf(paste0(workd,"/Results/Distance_Relationships_tumoradj_noTLS.pdf"),width=10,height=8)

##Normal
xx<-dist_Normal
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlistNormal<-colorlist[as.numeric(colnames(mat_Normal))] #set color
names(colorlistNormal)<-as.character(colnames(mat_Normal)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_Normal[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Normal",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_Normal,1/dist_CP,1/dist_PanIN,1/dist_PDAC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=0, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="lightgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistNormal[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype)), 
       col = cbbPalette , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

##CP
xx<-dist_CP
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlistCP<-colorlist[as.numeric(colnames(mat_CP))] #set color
names(colorlistCP)<-as.character(colnames(mat_CP)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_CP[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="CP",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_Normal,1/dist_CP,1/dist_PanIN,1/dist_PDAC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=0, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="lightgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistCP[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype)), 
       col = cbbPalette , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

##PanIN
xx<-dist_PanIN
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlistPanIN<-colorlist[as.numeric(colnames(mat_PanIN))] #set color
names(colorlistPanIN)<-as.character(colnames(mat_PanIN)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_PanIN[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="PanIN",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_Normal,1/dist_CP,1/dist_PanIN,1/dist_PDAC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=0, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="lightgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistPanIN[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype)), 
       col = cbbPalette , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

##PDAC
xx<-dist_PDAC
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlistPDAC<-colorlist[as.numeric(colnames(mat_PDAC))] #set color
names(colorlistPDAC)<-as.character(colnames(mat_PDAC)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_PDAC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="PDAC",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_Normal,1/dist_CP,1/dist_PanIN,1/dist_PDAC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=0, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="lightgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistPDAC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype)), 
       col = cbbPalette , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

dev.off()

