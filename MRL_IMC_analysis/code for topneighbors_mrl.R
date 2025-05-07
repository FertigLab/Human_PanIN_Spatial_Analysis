####NEIGHBORHOOD ANALYSIS####

## Set up dataframe for first and second closest objects to each index object
expr<-fsApply(output$fcs1, exprs)
expr<-data.frame(expr[,c("CellId","Num_Neighbor","NN1","NN2")],
                 objtype=output$cell_clustering1m,
                 sample_id=output$sample_ids)
rownames(expr)<-paste(expr$sample_id,expr$CellId,sep="_")

## Match the annotation to each neighbor object and annotate responses
# Comment out the samples not to be included in the analysis
# Samples not to be included are ones that are not biologically relevant and tumor samples that are not PDAC
expr$NN1 <- paste(expr$sample_id,expr$NN1,sep="_")
expr$NN2 <- paste(expr$sample_id,expr$NN2,sep="_")
expr$N1type<-expr[expr$NN1,]$objtype
expr$N2type<-expr[expr$NN2,]$objtype
expr$Type<- 1
expr$TLS_status<- 1
expr$Type[expr$sample_id %in% c("1C_28500_1_n_1",
                                "1C_28500_1_nad_1",
                                "1C_28500_2_n_1",
                                "1C_28500_2_nad_1",
                                "1C_28500_3_n_1",
                                "1C_28500_3_n_2",
                                "1C_28500_3_nad_1",
                                "1K_25386_1_n_1",
                                "1K_25386_1_nad_1",
                                "1K_25386_2_n_1",
                                "1K_25386_3_n_1",
                                "1K_25386_3_nad_1",
                                "1O_11554_5_n_1",
                                "1O_25590_5_n_1",
                                "2D_20447_6_n_1",
                                "2K_25590_2_n_1")]<- "Normal"
expr$Type[expr$sample_id %in% c(
                                # "1D_28500_4_pa_1",
                                "1D_28500_4_panin_1",
                                "1E_25386_5_panin_1",
                                "1H_25386_4_panin_1",
                                "1M_28500_1_panin_1",
                                "1M_28500_1_tls_panin_1",
                                "1M_28500_1_tls_panin_2",
                                # "1M_28500_2_pa_1",
                                "1M_28500_2_panin_1",
                                "1M_28500_2_panin_2",
                                # "1M_28500_3_pa_1",
                                "1M_28500_3_panin_1",
                                "1M_28500_3_tls_panin_1",
                                "1M_28500_3_tls_panin_2",
                                # "1M_28500_7_pa_1",
                                "1M_28500_7_panin_1",
                                "1O_11554_6_panin_1",
                                "1O_11554_6_tls_pa_1",
                                "1O_25590_6_panin_1",
                                "2D_20447_2_panin_1",
                                "2D_20447_3_panin_1",
                                "2D_20447_3_tls_panin_1",
                                # "2HH_20447_1_pa_1",
                                "2HH_20447_1_panin_1",
                                # "2O_20447_1_pa_1",
                                "2O_20447_1_panin_1",
                                # "2O_20447_2_pa_1",
                                # "2O_20447_2_pa_2",
                                "2O_20447_2_panin_1"
                                )]<-"PanIN"
expr$Type[expr$sample_id %in% c("1D_28500_3_cp_1",
                                "1E_25386_4_cp_1",
                                "1H_25386_1_cp_1",
                                "1O_11554_2_cp_1",
                                "1O_11554_2_tls_cp_1",
                                "1O_25590_2_cp_1",
                                "1O_25590_2_tls_cp_1",
                                "2D_20447_1_cp_1",
                                "2HH_20447_2_cp_1",
                                "2HH_20447_2_tls_cpb_1",
                                "2HH_20447_3_cp_1",
                                "2HH_20447_3_tls_cp_1",
                                "2HH_20447_3_tls_cp_2",
                                "2HH_20447_3_tls_cp_3",
                                "1M_28500_7_tls_panin_1"
                                )]<-"CP"
expr$Type[expr$sample_id %in% c(
                                # "1D_28500_1_tls_ta_1",
                                # "1D_28500_1_tls_ta_2",
                                # "1D_28500_1_tls_ta_3",
                                # "1D_28500_1_tls_ta_4",
                                # "1D_28500_2_t_1_1",
                                "1E_25386_1_t_1",
                                "1E_25386_2_ta_1",
                                "1E_25386_2_tls_ta_1",
                                "1E_25386_2_tls_ta_2",
                                "1E_25386_3_t_1",
                                "1H_25386_2_t_1",
                                # "1H_25386_3_tb_1",
                                # "1M_28500_5_t_1_1",
                                # "1M_28500_6_t_1_1",
                                "1O_11554_3_t_1",
                                "1O_11554_4_t_1",
                                # "1O_11554_4_tb_1",
                                "1O_25590_3_t_1",
                                "1O_25590_4_t_1",
                                # "1O_25590_4_tb_1",
                                "2D_20447_4_t_1",
                                # "2D_20447_4_tb_1",
                                "2D_20447_5_t_1",
                                "2K_25590_4_tls_tb_1",
                                # "2K_25590_5_tb_1",
                                "2O_20447_3_t_1",
                                "2O_20447_3_tls_tb_1",
                                "2O_20447_4_t_1",
                                "2O_20447_4_tls_t_1",
                                "2O_20447_4_tls_t_2",
                                "2O_20447_5_t_1",
                                "2O_20447_5_tls_t_1",
                                # "2O_20447_6_ta_1",
                                "2O_20447_6_tls_ta_1",
                                # "2O_20447_7_tb_1",
                                # "2O_20447_7_tb_2",
                                "2O_20447_7_tls_tb_1",
                                # "2O_20447_8_ta_1",
                                # "2O_20447_8_ta_2",
                                "2O_20447_8_tls_ta_1",
                                # "2O_20447_9_ta_1",
                                # "2O_20447_9_tb_1",
                                "2O_20447_9_tls_ta_1",
                                # "2O_20447_10_tb_1",
                                "2O_20447_10_tls_tb_1",
                                # "2Q_25590_1_tb_1",
                                # "2Q_25590_2_ta_1",
                                "2Q_25590_2_tls_ta_1",
                                "2Q_25590_2_tls_ta_2",
                                "2Q_25590_2_tls_ta_3"
                                )]<-"PDAC"
expr$TLS_status[expr$sample_id %in% c(
                                      # "1D_28500_1_tls_ta_1",
                                      # "1D_28500_1_tls_ta_2",
                                      # "1D_28500_1_tls_ta_3",
                                      # "1D_28500_1_tls_ta_4",
                                      "1E_25386_2_tls_ta_1",
                                      "1E_25386_2_tls_ta_2",
                                      "1M_28500_1_tls_panin_1",
                                      "1M_28500_1_tls_panin_2",
                                      "1M_28500_3_tls_panin_1",
                                      "1M_28500_3_tls_panin_2",
                                      "1M_28500_7_tls_panin_1",
                                      "1O_11554_2_tls_cp_1",
                                      "1O_11554_6_tls_pa_1",
                                      "1O_25590_2_tls_cp_1",
                                      "2D_20447_3_tls_panin_1",
                                      "2HH_20447_2_tls_cpb_1",
                                      "2HH_20447_3_tls_cp_1",
                                      "2HH_20447_3_tls_cp_2",
                                      "2HH_20447_3_tls_cp_3",
                                      "2K_25590_1_tls_1",
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
                                      "2Q_25590_2_tls_ta_3")]<-"TLS"
expr$TLS_status[expr$sample_id %in% c("1C_28500_1_n_1",
                                      "1C_28500_1_nad_1",
                                      "1C_28500_2_n_1",
                                      "1C_28500_2_nad_1",
                                      "1C_28500_3_n_1",
                                      "1C_28500_3_n_2",
                                      "1C_28500_3_nad_1",
                                      "1D_28500_2_t_1_1",
                                      "1D_28500_3_cp_1",
                                      # "1D_28500_4_pa_1",
                                      "1D_28500_4_panin_1",
                                      "1E_25386_1_t_1",
                                      # "1E_25386_2_ta_1",
                                      "1E_25386_3_t_1",
                                      "1E_25386_4_cp_1",
                                      "1E_25386_5_panin_1",
                                      "1H_25386_1_cp_1",
                                      "1H_25386_2_t_1",
                                      # "1H_25386_3_tb_1",
                                      "1H_25386_4_panin_1",
                                      "1K_25386_1_n_1",
                                      "1K_25386_1_nad_1",
                                      "1K_25386_2_n_1",
                                      "1K_25386_3_n_1",
                                      "1K_25386_3_nad_1",
                                      "1M_28500_1_panin_1",
                                      # "1M_28500_2_pa_1",
                                      "1M_28500_2_panin_1",
                                      "1M_28500_2_panin_2",
                                      # "1M_28500_3_pa_1",
                                      "1M_28500_3_panin_1",
                                      # "1M_28500_5_t_1_1",
                                      # "1M_28500_6_t_1_1",
                                      # "1M_28500_7_pa_1",
                                      "1M_28500_7_panin_1",
                                      "1O_11554_2_cp_1",
                                      "1O_11554_3_t_1",
                                      "1O_11554_4_t_1",
                                      # "1O_11554_4_tb_1",
                                      "1O_11554_5_n_1",
                                      "1O_11554_6_panin_1",
                                      "1O_25590_2_cp_1",
                                      "1O_25590_3_t_1",
                                      "1O_25590_4_t_1",
                                      # "1O_25590_4_tb_1",
                                      "1O_25590_5_n_1",
                                      "1O_25590_6_panin_1",
                                      "2D_20447_1_cp_1",
                                      "2D_20447_2_panin_1",
                                      "2D_20447_3_panin_1",
                                      "2D_20447_4_t_1",
                                      # "2D_20447_4_tb_1",
                                      "2D_20447_5_t_1",
                                      "2D_20447_6_n_1",
                                      # "2HH_20447_1_pa_1",
                                      "2HH_20447_1_panin_1",
                                      "2HH_20447_2_cp_1",
                                      "2HH_20447_3_cp_1",
                                      "2K_25590_2_n_1",
                                      # "2K_25590_5_tb_1",
                                      # "2O_20447_1_pa_1",
                                      "2O_20447_1_panin_1",
                                      # "2O_20447_2_pa_1",
                                      # "2O_20447_2_pa_2",
                                      "2O_20447_2_panin_1",
                                      "2O_20447_3_t_1",
                                      "2O_20447_4_t_1",
                                      "2O_20447_5_t_1"
                                      # "2O_20447_6_ta_1",
                                      # "2O_20447_7_tb_1",
                                      # "2O_20447_7_tb_2",
                                      # "2O_20447_8_ta_1",
                                      # "2O_20447_8_ta_2",
                                      # "2O_20447_9_ta_1",
                                      # "2O_20447_9_tb_1",
                                      # "2O_20447_10_tb_1",
                                      # "2Q_25590_1_tb_1",
                                      # "2Q_25590_2_ta_1"
                                      )]<-"nonTLS"




# First subset the dataframe to removed unassigned cells (cluster annotation is "Unassigned") as well as nonlymphoid cell subtypes
expr_assigned <-expr[expr$objtype!="Unassigned",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CD68+ CD16+ CD11c+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CD68+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="DCSIGN+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CD68+ CD16+ HLADR+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CK+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="PDPN+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="COL+ SMA+ VIM+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CD57+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CK+ PDL1+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CD57+ TOX2+",]
expr_assigned <-expr_assigned[expr_assigned$objtype!="CD45+",]


# Remove cells with nearest neighbors that are unassigned cells as well as nonlymphoid cell subtypes
expr_assigned1 <- expr_assigned[expr_assigned$N1type!="Unassigned",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CD68+ CD16+ CD11c+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CD68+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="DCSIGN+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CD68+ CD16+ HLADR+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CK+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="PDPN+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="COL+ SMA+ VIM+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CD57+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CK+ PDL1+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CD57+ TOX2+",]
expr_assigned1 <-expr_assigned1[expr_assigned1$N1type!="CD45+",]

expr_assigned2 <- expr_assigned1[expr_assigned1$N2type!="Unassigned",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CD68+ CD16+ CD11c+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CD68+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="DCSIGN+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CD68+ CD16+ HLADR+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CK+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="PDPN+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="COL+ SMA+ VIM+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CD57+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CK+ PDL1+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CD57+ TOX2+",]
expr_assigned2 <-expr_assigned2[expr_assigned2$N2type!="CD45+",]


# If you don't want to subset any cells out of the analysis
expr_assigned2 <- expr

## Make a heatmap of cumulative neighbor data (top 2 neighbors for each object)
nb_expr<-expr_assigned2[expr_assigned2$Num_Neighbor!=0,]
Heatmap(as.matrix.data.frame(table(nb_expr$objtype,nb_expr$N1type)))
df_nb_expr<-data.frame(objtype=nb_expr$objtype,N1type=nb_expr$N1type)
m_nb_expr<-as.matrix.data.frame(table(df_nb_expr$objtype,df_nb_expr$N1type), rownames.force = T)
Heatmap(m_nb_expr)

## Subset dataframes by type of pathology
nb_expr_n<-as.data.frame(nb_expr[nb_expr$Type=="Normal",])
nb_expr_panin<-as.data.frame(nb_expr[nb_expr$Type=="PanIN",])
nb_expr_cp<-as.data.frame(nb_expr[nb_expr$Type=="CP",])
nb_expr_pdac<-as.data.frame(nb_expr[nb_expr$Type=="PDAC",])

## Subset these dataframes by another level to separate TLS vs nonTLS samples in each region

nb_expr_n_tls<-as.data.frame(nb_expr_n[nb_expr_n$TLS_status=="TLS",])
nb_expr_n_nontls<-as.data.frame(nb_expr_n[nb_expr_n$TLS_status=="nonTLS",])

nb_expr_panin_tls<-as.data.frame(nb_expr_panin[nb_expr_panin$TLS_status=="TLS",])
nb_expr_panin_nontls<-as.data.frame(nb_expr_panin[nb_expr_panin$TLS_status=="nonTLS",])

nb_expr_cp_tls<-as.data.frame(nb_expr_cp[nb_expr_cp$TLS_status=="TLS",])
nb_expr_cp_nontls<-as.data.frame(nb_expr_cp[nb_expr_cp$TLS_status=="nonTLS",])

nb_expr_pdac_tls<-as.data.frame(nb_expr_pdac[nb_expr_pdac$TLS_status=="TLS",])
nb_expr_pdac_nontls<-as.data.frame(nb_expr_pdac[nb_expr_pdac$TLS_status=="nonTLS",])

## Make heatmaps for each type and TLS status for the first closest object
# First for normal samples, combined data as well as dataframes separating TLS and nonTLS (in this case there are no TLSs in normal samples)
m_nb_expr_n<-as.matrix.data.frame(table(nb_expr_n$objtype,nb_expr_n$N1type), rownames.force = T)
colnames(m_nb_expr_n)<-rownames(m_nb_expr_n)
Heatmap(m_nb_expr_n,cluster_rows = F,cluster_columns = F)

m_nb_expr_n_tls<-as.matrix.data.frame(table(nb_expr_n_tls$objtype,nb_expr_n_tls$N1type), rownames.force = T)
colnames(m_nb_expr_n_tls)<-rownames(m_nb_expr_n_tls)
Heatmap(m_nb_expr_n_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_n_nontls<-as.matrix.data.frame(table(nb_expr_n_nontls$objtype,nb_expr_n_nontls$N1type), rownames.force = T)
colnames(m_nb_expr_n_nontls)<-rownames(m_nb_expr_n_nontls)
Heatmap(m_nb_expr_n_nontls,cluster_rows = F,cluster_columns = F)

# For PanIN samples, combined data as well as dataframes separating TLS and nonTLS 
m_nb_expr_panin<-as.matrix.data.frame(table(nb_expr_panin$objtype,nb_expr_panin$N1type), rownames.force = T)
colnames(m_nb_expr_panin)<-rownames(m_nb_expr_panin)
Heatmap(m_nb_expr_panin,cluster_rows = F,cluster_columns = F)

m_nb_expr_panin_tls<-as.matrix.data.frame(table(nb_expr_panin_tls$objtype,nb_expr_panin_tls$N1type), rownames.force = T)
colnames(m_nb_expr_panin_tls)<-rownames(m_nb_expr_panin_tls)
Heatmap(m_nb_expr_panin_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_panin_nontls<-as.matrix.data.frame(table(nb_expr_panin_nontls$objtype,nb_expr_panin_nontls$N1type), rownames.force = T)
colnames(m_nb_expr_panin_nontls)<-rownames(m_nb_expr_panin_nontls)
Heatmap(m_nb_expr_panin_nontls,cluster_rows = F,cluster_columns = F)

# For CP samples, combined data as well as dataframes separating TLS and nonTLS 
m_nb_expr_cp<-as.matrix.data.frame(table(nb_expr_cp$objtype,nb_expr_cp$N1type), rownames.force = T)
colnames(m_nb_expr_cp)<-rownames(m_nb_expr_cp)
Heatmap(m_nb_expr_cp,cluster_rows = F,cluster_columns = F)

m_nb_expr_cp_tls<-as.matrix.data.frame(table(nb_expr_cp_tls$objtype,nb_expr_cp_tls$N1type), rownames.force = T)
colnames(m_nb_expr_cp_tls)<-rownames(m_nb_expr_cp_tls)
Heatmap(m_nb_expr_cp_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_cp_nontls<-as.matrix.data.frame(table(nb_expr_cp_nontls$objtype,nb_expr_cp_nontls$N1type), rownames.force = T)
colnames(m_nb_expr_cp_nontls)<-rownames(m_nb_expr_cp_nontls)
Heatmap(m_nb_expr_cp_nontls,cluster_rows = F,cluster_columns = F)

# For PDAC samples, combined data as well as dataframes separating TLS and nonTLS 
m_nb_expr_pdac<-as.matrix.data.frame(table(nb_expr_pdac$objtype,nb_expr_pdac$N1type), rownames.force = T)
colnames(m_nb_expr_pdac)<-rownames(m_nb_expr_pdac)
Heatmap(m_nb_expr_pdac,cluster_rows = F,cluster_columns = F)

m_nb_expr_pdac_tls<-as.matrix.data.frame(table(nb_expr_pdac_tls$objtype,nb_expr_pdac_tls$N1type), rownames.force = T)
colnames(m_nb_expr_pdac_tls)<-rownames(m_nb_expr_pdac_tls)
Heatmap(m_nb_expr_pdac_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_pdac_nontls<-as.matrix.data.frame(table(nb_expr_pdac_nontls$objtype,nb_expr_pdac_nontls$N1type), rownames.force = T)
colnames(m_nb_expr_pdac_nontls)<-rownames(m_nb_expr_pdac_nontls)
Heatmap(m_nb_expr_pdac_nontls,cluster_rows = F,cluster_columns = F)

## Make heatmaps for each sample based on type and tls status for the second closest object
# For combined normal samples, and then separated TLS and nonTLS samples
m_nb_expr_n2<-as.matrix.data.frame(table(nb_expr_n$objtype,nb_expr_n$N2type), rownames.force = T)
colnames(m_nb_expr_n2)<-rownames(m_nb_expr_n2)
Heatmap(m_nb_expr_n2,cluster_rows = F,cluster_columns = F)

m_nb_expr_n2_tls<-as.matrix.data.frame(table(nb_expr_n_tls$objtype,nb_expr_n_tls$N2type), rownames.force = T)
colnames(m_nb_expr_n2_tls)<-rownames(m_nb_expr_n2_tls)
Heatmap(m_nb_expr_n2_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_n2_nontls<-as.matrix.data.frame(table(nb_expr_n_nontls$objtype,nb_expr_n_nontls$N2type), rownames.force = T)
colnames(m_nb_expr_n2_nontls)<-rownames(m_nb_expr_n2_nontls)
Heatmap(m_nb_expr_n2,cluster_rows = F,cluster_columns = F)

# For combined PanIN samples, and then separated TLS and nonTLS samples
m_nb_expr_panin2<-as.matrix.data.frame(table(nb_expr_panin$objtype,nb_expr_panin$N2type), rownames.force = T)
colnames(m_nb_expr_panin2)<-rownames(m_nb_expr_panin2)
Heatmap(m_nb_expr_panin2,cluster_rows = F,cluster_columns = F)

m_nb_expr_panin2_tls<-as.matrix.data.frame(table(nb_expr_panin_tls$objtype,nb_expr_panin_tls$N2type), rownames.force = T)
colnames(m_nb_expr_panin2_tls)<-rownames(m_nb_expr_panin2_tls)
Heatmap(m_nb_expr_panin2_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_panin2_nontls<-as.matrix.data.frame(table(nb_expr_panin_nontls$objtype,nb_expr_panin_nontls$N2type), rownames.force = T)
colnames(m_nb_expr_panin2_nontls)<-rownames(m_nb_expr_panin2_nontls)
Heatmap(m_nb_expr_panin2_nontls,cluster_rows = F,cluster_columns = F)

# For combined CP samples, and then separated TLS and nonTLS samples
m_nb_expr_cp2<-as.matrix.data.frame(table(nb_expr_cp$objtype,nb_expr_cp$N2type), rownames.force = T)
colnames(m_nb_expr_cp2)<-rownames(m_nb_expr_cp2)
Heatmap(m_nb_expr_cp2,cluster_rows = F,cluster_columns = F)

m_nb_expr_cp2_tls<-as.matrix.data.frame(table(nb_expr_cp_tls$objtype,nb_expr_cp_tls$N2type), rownames.force = T)
colnames(m_nb_expr_cp2_tls)<-rownames(m_nb_expr_cp2_tls)
Heatmap(m_nb_expr_cp2_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_cp2_nontls<-as.matrix.data.frame(table(nb_expr_cp_nontls$objtype,nb_expr_cp_nontls$N2type), rownames.force = T)
colnames(m_nb_expr_cp2_nontls)<-rownames(m_nb_expr_cp2_nontls)
Heatmap(m_nb_expr_cp2_nontls,cluster_rows = F,cluster_columns = F)

# For combined PDAC samples, and then separated TLS and nonTLS samples
m_nb_expr_pdac2<-as.matrix.data.frame(table(nb_expr_pdac$objtype,nb_expr_pdac$N2type), rownames.force = T)
colnames(m_nb_expr_pdac2)<-rownames(m_nb_expr_pdac2)
Heatmap(m_nb_expr_pdac2,cluster_rows = F,cluster_columns = F)

m_nb_expr_pdac2_tls<-as.matrix.data.frame(table(nb_expr_pdac_tls$objtype,nb_expr_pdac_tls$N2type), rownames.force = T)
colnames(m_nb_expr_pdac2_tls)<-rownames(m_nb_expr_pdac2_tls)
Heatmap(m_nb_expr_pdac2_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_pdac2_nontls<-as.matrix.data.frame(table(nb_expr_pdac_nontls$objtype,nb_expr_pdac_nontls$N2type), rownames.force = T)
colnames(m_nb_expr_pdac2_nontls)<-rownames(m_nb_expr_pdac2_nontls)
Heatmap(m_nb_expr_pdac2_nontls,cluster_rows = F,cluster_columns = F)

## Make heatmaps for all tissue types and tls status for both first and second closest object
# Combined normal, tls and non tls first and second neighbor
m_nb_expr_ntop2<-as.matrix.data.frame(table(nb_expr_n$objtype,nb_expr_n$N1type)+table(nb_expr_n$objtype,nb_expr_n$N2type), rownames.force = T)
colnames(m_nb_expr_ntop2)<-rownames(m_nb_expr_ntop2)
Heatmap(m_nb_expr_ntop2,cluster_rows = F,cluster_columns = F)

m_nb_expr_ntop2_tls<-as.matrix.data.frame(table(nb_expr_n_tls$objtype,nb_expr_n_tls$N1type)+table(nb_expr_n_tls$objtype,nb_expr_n_tls$N2type), rownames.force = T)
colnames(m_nb_expr_ntop2_tls)<-rownames(m_nb_expr_ntop2_tls)
Heatmap(m_nb_expr_ntop2_tls,cluster_rows = F,cluster_columns = F)

m_nb_expr_ntop2_nontls<-as.matrix.data.frame(table(nb_expr_n_nontls$objtype,nb_expr_n_nontls$N1type)+table(nb_expr_n_nontls$objtype,nb_expr_n_nontls$N2type), rownames.force = T)
m_nb_expr_ntop2_nontls2 <- m_nb_expr_ntop2_nontls[c("CD20+ CD21+ CD23+",
                                                    "CD20+ CD45RA+",
                                                    "CD4+ CD45RO+ CCR7+ HLADR+",
                                                    "CD8+ GZMB+",
                                                    "CD4+ CD45RO+",
                                                    "CD8+ CD45RO+",
                                                    "CD4+ CD45RO+ TOX2+ PD1+",
                                                    "CD4+ CD45RA+",
                                                    "CD4+ CD45RA+ KI67+",
                                                    "CD8+",
                                                    "CD4+ FOXP3+ PD1-",
                                                    "CD8+ GZMB+ KI67+ LAG3+"
                                                    # "CD68+ CD16+ CD11c+",
                                                    # "CD68+",
                                                    # "DCSIGN+",
                                                    # "CD68+ CD16+ HLADR+"
                                                    # "CK+",
                                                    # "CK+ PDL1+"
                                                    ),]

colnames(m_nb_expr_ntop2_nontls2)<-rownames(m_nb_expr_ntop2_nontls)

pdf('m_nb_expr_ntop2_nontls2_pdacpts.pdf',width=5,height=3.5)

Heatmap(m_nb_expr_ntop2_nontls2[,c("CD20+ CD21+ CD23+",
                                   "CD20+ CD45RA+",
                                   "CD4+ CD45RO+ CCR7+ HLADR+",
                                   "CD8+ GZMB+",
                                   "CD4+ CD45RO+",
                                   "CD8+ CD45RO+",
                                   "CD4+ CD45RO+ TOX2+ PD1+",
                                   "CD4+ CD45RA+",
                                   "CD4+ CD45RA+ KI67+",
                                   "CD8+",
                                   "CD4+ FOXP3+ PD1-",
                                   "CD8+ GZMB+ KI67+ LAG3+"
                                   # "CD68+ CD16+ CD11c+",
                                   # "CD68+",
                                   # "DCSIGN+",
                                   # "CD68+ CD16+ HLADR+"
                                   # "CK+",
                                   # "CK+ PDL1+"
                                   )], 
                                    row_names_gp = gpar(fontsize = 6),
                                    column_names_gp = gpar(fontsize = 6),
                                    cluster_rows = T,
                                    cluster_columns = T,
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    heatmap_legend_param = list(title = "Cells"))
dev.off()

# Combined PanIN, tls and non tls first and second neighbor
m_nb_expr_panintop2<-as.matrix.data.frame(table(nb_expr_panin$objtype,nb_expr_panin$N1type)+table(nb_expr_panin$objtype,nb_expr_panin$N2type), rownames.force = T)
colnames(m_nb_expr_panintop2)<-rownames(m_nb_expr_panintop2)
Heatmap(m_nb_expr_panintop2,cluster_rows = F,cluster_columns = F)

m_nb_expr_panintop2_tls<-as.matrix.data.frame(table(nb_expr_panin_tls$objtype,nb_expr_panin_tls$N1type)+table(nb_expr_panin_tls$objtype,nb_expr_panin_tls$N2type), rownames.force = T)
m_nb_expr_panintop2_tls2 <- m_nb_expr_panintop2_tls[c("CD20+ CD21+ CD23+",
                                                      "CD20+ CD45RA+",
                                                      "CD4+ CD45RO+ CCR7+ HLADR+",
                                                      "CD8+ GZMB+",
                                                      "CD4+ CD45RO+",
                                                      "CD8+ CD45RO+",
                                                      "CD4+ CD45RO+ TOX2+ PD1+",
                                                      "CD4+ CD45RA+",
                                                      "CD4+ CD45RA+ KI67+",
                                                      "CD8+",
                                                      "CD4+ FOXP3+ PD1-",
                                                      "CD8+ GZMB+ KI67+ LAG3+"
                                                      # "CD68+ CD16+ CD11c+",
                                                      # "CD68+",
                                                      # "DCSIGN+",
                                                      # "CD68+ CD16+ HLADR+"
                                                      ),]

colnames(m_nb_expr_panintop2_tls2)<-rownames(m_nb_expr_panintop2_tls)

pdf('m_nb_expr_panintop2_tls2_pdacpts.pdf',width=5,height=3.5)

Heatmap(m_nb_expr_panintop2_tls2[,c("CD20+ CD21+ CD23+",
                                    "CD20+ CD45RA+",
                                    "CD4+ CD45RO+ CCR7+ HLADR+",
                                    "CD8+ GZMB+",
                                    "CD4+ CD45RO+",
                                    "CD8+ CD45RO+",
                                    "CD4+ CD45RO+ TOX2+ PD1+",
                                    "CD4+ CD45RA+",
                                    "CD4+ CD45RA+ KI67+",
                                    "CD8+",
                                    "CD4+ FOXP3+ PD1-",
                                    "CD8+ GZMB+ KI67+ LAG3+"
                                    # "CD68+ CD16+ CD11c+",
                                    # "CD68+",
                                    # "DCSIGN+",
                                    # "CD68+ CD16+ HLADR+"
                                    )], 
                                    row_names_gp = gpar(fontsize = 6),
                                    column_names_gp = gpar(fontsize = 6),
                                    cluster_rows = T,
                                    cluster_columns = T,
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    heatmap_legend_param = list(title = "Cells"))

dev.off()

m_nb_expr_panintop2_nontls<-as.matrix.data.frame(table(nb_expr_panin_nontls$objtype,nb_expr_panin_nontls$N1type)+table(nb_expr_panin_nontls$objtype,nb_expr_panin_nontls$N2type), rownames.force = T)
m_nb_expr_panintop2_nontls2 <- m_nb_expr_panintop2_nontls[c("CD20+ CD21+ CD23+",
                                                            "CD20+ CD45RA+",
                                                            "CD4+ CD45RO+ CCR7+ HLADR+",
                                                            "CD8+ GZMB+",
                                                            "CD4+ CD45RO+",
                                                            "CD8+ CD45RO+",
                                                            "CD4+ CD45RO+ TOX2+ PD1+",
                                                            "CD4+ CD45RA+",
                                                            "CD4+ CD45RA+ KI67+",
                                                            "CD8+",
                                                            "CD4+ FOXP3+ PD1-",
                                                            "CD8+ GZMB+ KI67+ LAG3+"
                                                            # "CD68+ CD16+ CD11c+",
                                                            # "CD68+",
                                                            # "DCSIGN+",
                                                            # "CD68+ CD16+ HLADR+"
                                                            ),]
colnames(m_nb_expr_panintop2_nontls2)<-rownames(m_nb_expr_panintop2_nontls)

pdf('m_nb_expr_panintop2_nontls2_pdacpts.pdf',width=5,height=3.5)

Heatmap(m_nb_expr_panintop2_nontls2[,c("CD20+ CD21+ CD23+",
                                       "CD20+ CD45RA+",
                                       "CD4+ CD45RO+ CCR7+ HLADR+",
                                       "CD8+ GZMB+",
                                       "CD4+ CD45RO+",
                                       "CD8+ CD45RO+",
                                       "CD4+ CD45RO+ TOX2+ PD1+",
                                       "CD4+ CD45RA+",
                                       "CD4+ CD45RA+ KI67+",
                                       "CD8+",
                                       "CD4+ FOXP3+ PD1-",
                                       "CD8+ GZMB+ KI67+ LAG3+"
                                       # "CD68+ CD16+ CD11c+",
                                       # "CD68+",
                                       # "DCSIGN+",
                                       # "CD68+ CD16+ HLADR+"
                                       )], 
                                        row_names_gp = gpar(fontsize = 6),
                                        column_names_gp = gpar(fontsize = 6),
                                        cluster_rows = T,
                                        cluster_columns = T,
                                        show_row_dend = F,
                                        show_column_dend = F,
                                        heatmap_legend_param = list(title = "Cells"))

dev.off()

# Combined CP, tls and non tls first and second neighbor
m_nb_expr_cptop2<-as.matrix.data.frame(table(nb_expr_cp$objtype,nb_expr_cp$N1type)+table(nb_expr_cp$objtype,nb_expr_cp$N2type), rownames.force = T)
colnames(m_nb_expr_cptop2)<-rownames(m_nb_expr_cptop2)
Heatmap(m_nb_expr_cptop2,cluster_rows = F,cluster_columns = F)

m_nb_expr_cptop2_tls<-as.matrix.data.frame(table(nb_expr_cp_tls$objtype,nb_expr_cp_tls$N1type)+table(nb_expr_cp_tls$objtype,nb_expr_cp_tls$N2type), rownames.force = T)
m_nb_expr_cptop2_tls2 <- m_nb_expr_cptop2_tls[c("CD20+ CD21+ CD23+",
                                                "CD20+ CD45RA+",
                                                "CD4+ CD45RO+ CCR7+ HLADR+",
                                                "CD8+ GZMB+",
                                                "CD4+ CD45RO+",
                                                "CD8+ CD45RO+",
                                                "CD4+ CD45RO+ TOX2+ PD1+",
                                                "CD4+ CD45RA+",
                                                "CD4+ CD45RA+ KI67+",
                                                "CD8+",
                                                "CD4+ FOXP3+ PD1-",
                                                "CD8+ GZMB+ KI67+ LAG3+"
                                                # "CD68+ CD16+ CD11c+",
                                                # "CD68+",
                                                # "DCSIGN+",
                                                # "CD68+ CD16+ HLADR+"
                                                ),]
colnames(m_nb_expr_cptop2_tls2)<-rownames(m_nb_expr_cptop2_tls)

pdf('m_nb_expr_cptop2_tls2_pdacpts.pdf',width=5,height=3.5)

Heatmap(m_nb_expr_cptop2_tls2[,c("CD20+ CD21+ CD23+",
                                 "CD20+ CD45RA+",
                                 "CD4+ CD45RO+ CCR7+ HLADR+",
                                 "CD8+ GZMB+",
                                 "CD4+ CD45RO+",
                                 "CD8+ CD45RO+",
                                 "CD4+ CD45RO+ TOX2+ PD1+",
                                 "CD4+ CD45RA+",
                                 "CD4+ CD45RA+ KI67+",
                                 "CD8+",
                                 "CD4+ FOXP3+ PD1-",
                                 "CD8+ GZMB+ KI67+ LAG3+"
                                 # "CD68+ CD16+ CD11c+",
                                 # "CD68+",
                                 # "DCSIGN+",
                                 # "CD68+ CD16+ HLADR+"
                                 )], 
                                  row_names_gp = gpar(fontsize = 6),
                                  column_names_gp = gpar(fontsize = 6),
                                  cluster_rows = T,
                                  cluster_columns = T,
                                  show_row_dend = F,
                                  show_column_dend = F,
                                  heatmap_legend_param = list(title = "Cells"))

dev.off()

m_nb_expr_cptop2_nontls<-as.matrix.data.frame(table(nb_expr_cp_nontls$objtype,nb_expr_cp_nontls$N1type)+table(nb_expr_cp_nontls$objtype,nb_expr_cp_nontls$N2type), rownames.force = T)
m_nb_expr_cptop2_nontls2 <- m_nb_expr_cptop2_nontls[c("CD20+ CD21+ CD23+",
                                                      "CD20+ CD45RA+",
                                                      "CD4+ CD45RO+ CCR7+ HLADR+",
                                                      "CD8+ GZMB+",
                                                      "CD4+ CD45RO+",
                                                      "CD8+ CD45RO+",
                                                      "CD4+ CD45RO+ TOX2+ PD1+",
                                                      "CD4+ CD45RA+",
                                                      "CD4+ CD45RA+ KI67+",
                                                      "CD8+",
                                                      "CD4+ FOXP3+ PD1-",
                                                      "CD8+ GZMB+ KI67+ LAG3+"
                                                      # "CD68+ CD16+ CD11c+",
                                                      # "CD68+",
                                                      # "DCSIGN+",
                                                      # "CD68+ CD16+ HLADR+"
                                                      ),]
colnames(m_nb_expr_cptop2_nontls2)<-rownames(m_nb_expr_cptop2_nontls)

pdf('m_nb_expr_cptop2_nontls2_pdacpts.pdf',width=5,height=3.5)

Heatmap(m_nb_expr_cptop2_nontls2[,c("CD20+ CD21+ CD23+",
                                    "CD20+ CD45RA+",
                                    "CD4+ CD45RO+ CCR7+ HLADR+",
                                    "CD8+ GZMB+",
                                    "CD4+ CD45RO+",
                                    "CD8+ CD45RO+",
                                    "CD4+ CD45RO+ TOX2+ PD1+",
                                    "CD4+ CD45RA+",
                                    "CD4+ CD45RA+ KI67+",
                                    "CD8+",
                                    "CD4+ FOXP3+ PD1-",
                                    "CD8+ GZMB+ KI67+ LAG3+"
                                    # "CD68+ CD16+ CD11c+",
                                    # "CD68+",
                                    # "DCSIGN+",
                                    # "CD68+ CD16+ HLADR+"
                                    )], 
                                    row_names_gp = gpar(fontsize = 6),
                                    column_names_gp = gpar(fontsize = 6),
                                    cluster_rows = T,
                                    cluster_columns = T,
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    heatmap_legend_param = list(title = "Cells"))

dev.off()

# Combined PDAC, tls and non tls first and second neighbor
m_nb_expr_pdactop2<-as.matrix.data.frame(table(nb_expr_pdac$objtype,nb_expr_pdac$N1type)+table(nb_expr_pdac$objtype,nb_expr_pdac$N2type), rownames.force = T)
colnames(m_nb_expr_pdactop2)<-rownames(m_nb_expr_pdactop2)
Heatmap(m_nb_expr_pdactop2,cluster_rows = F,cluster_columns = F)

m_nb_expr_pdactop2_tls<-as.matrix.data.frame(table(nb_expr_pdac_tls$objtype,nb_expr_pdac_tls$N1type)+table(nb_expr_pdac_tls$objtype,nb_expr_pdac_tls$N2type), rownames.force = T)
m_nb_expr_pdactop2_tls2 <- m_nb_expr_pdactop2_tls[c("CD20+ CD21+ CD23+",
                                                    "CD20+ CD45RA+",
                                                    "CD4+ CD45RO+ CCR7+ HLADR+",
                                                    "CD8+ GZMB+",
                                                    "CD4+ CD45RO+",
                                                    "CD8+ CD45RO+",
                                                    "CD4+ CD45RO+ TOX2+ PD1+",
                                                    "CD4+ CD45RA+",
                                                    "CD4+ CD45RA+ KI67+",
                                                    "CD8+",
                                                    "CD4+ FOXP3+ PD1-",
                                                    "CD8+ GZMB+ KI67+ LAG3+"
                                                    # "CD68+ CD16+ CD11c+",
                                                    # "CD68+",
                                                    # "DCSIGN+",
                                                    # "CD68+ CD16+ HLADR+"
                                                    ),]
colnames(m_nb_expr_pdactop2_tls2)<-rownames(m_nb_expr_pdactop2_tls)

pdf('m_nb_expr_pdactop2_tls2_pdacpts.pdf',width=5,height=3.5)

Heatmap(m_nb_expr_pdactop2_tls2[,c("CD20+ CD21+ CD23+",
                                   "CD20+ CD45RA+",
                                   "CD4+ CD45RO+ CCR7+ HLADR+",
                                   "CD8+ GZMB+",
                                   "CD4+ CD45RO+",
                                   "CD8+ CD45RO+",
                                   "CD4+ CD45RO+ TOX2+ PD1+",
                                   "CD4+ CD45RA+",
                                   "CD4+ CD45RA+ KI67+",
                                   "CD8+",
                                   "CD4+ FOXP3+ PD1-",
                                   "CD8+ GZMB+ KI67+ LAG3+"
                                   # "CD68+ CD16+ CD11c+",
                                   # "CD68+",
                                   # "DCSIGN+",
                                   # "CD68+ CD16+ HLADR+"
                                   )], 
                                    row_names_gp = gpar(fontsize = 6),
                                    column_names_gp = gpar(fontsize = 6),
                                    cluster_rows = T,
                                    cluster_columns = T,
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    heatmap_legend_param = list(title = "Cells"))

dev.off()

m_nb_expr_pdactop2_nontls<-as.matrix.data.frame(table(nb_expr_pdac_nontls$objtype,nb_expr_pdac_nontls$N1type)+table(nb_expr_pdac_nontls$objtype,nb_expr_pdac_nontls$N2type), rownames.force = T)
m_nb_expr_pdactop2_nontls2 <- m_nb_expr_pdactop2_nontls[c("CD20+ CD21+ CD23+",
                                                          "CD20+ CD45RA+",
                                                          "CD4+ CD45RO+ CCR7+ HLADR+",
                                                          "CD8+ GZMB+",
                                                          "CD4+ CD45RO+",
                                                          "CD8+ CD45RO+",
                                                          "CD4+ CD45RO+ TOX2+ PD1+",
                                                          "CD4+ CD45RA+",
                                                          "CD4+ CD45RA+ KI67+",
                                                          "CD8+",
                                                          "CD4+ FOXP3+ PD1-",
                                                          "CD8+ GZMB+ KI67+ LAG3+"
                                                          # "CD68+ CD16+ CD11c+",
                                                          # "CD68+",
                                                          # "DCSIGN+",
                                                          # "CD68+ CD16+ HLADR+"
                                                          ),]
colnames(m_nb_expr_pdactop2_nontls2)<-rownames(m_nb_expr_pdactop2_nontls)

pdf('m_nb_expr_pdactop2_nontls2_pdacpts.pdf',width=5,height=3.5)

Heatmap(m_nb_expr_pdactop2_nontls2[,c("CD20+ CD21+ CD23+",
                                      "CD20+ CD45RA+",
                                      "CD4+ CD45RO+ CCR7+ HLADR+",
                                      "CD8+ GZMB+",
                                      "CD4+ CD45RO+",
                                      "CD8+ CD45RO+",
                                      "CD4+ CD45RO+ TOX2+ PD1+",
                                      "CD4+ CD45RA+",
                                      "CD4+ CD45RA+ KI67+",
                                      "CD8+",
                                      "CD4+ FOXP3+ PD1-",
                                      "CD8+ GZMB+ KI67+ LAG3+"
                                      # "CD68+ CD16+ CD11c+",
                                      # "CD68+",
                                      # "DCSIGN+",
                                      # "CD68+ CD16+ HLADR+"
                                      )], 
                                      row_names_gp = gpar(fontsize = 6),
                                      column_names_gp = gpar(fontsize = 6),
                                      cluster_rows = T,
                                      cluster_columns = T,
                                      show_row_dend = F,
                                      show_column_dend = F,
                                      heatmap_legend_param = list(title = "Cells"))

dev.off()

## Make a differential heatmap (PDAC - PanIN)
Heatmap(m_nb_expr_pdactop2-m_nb_expr_panintop2, cluster_rows = F, cluster_columns = F)

## Annotations for optional labeling
ra<-rowAnnotation(clusters=names(colorassigned), col=list(clusters=colorassigned))
ca<-columnAnnotation(clusters=names(colorassigned), col=list(clusters=colorassigned))

## Heatmap of the differential heatmap
hp <- Heatmap(m_nb_expr_pdactop2-m_nb_expr_panintop2, 
              cluster_rows = F, cluster_columns = F,
              row_names_side="left",
              row_dend_side="right",
              row_title = "Index Clusters",
              column_title = "Top Neighboring Clusters",
              column_names_rot = 45,
              row_names_gp=gpar(fontsize=7),
              column_names_gp=gpar(fontsize=7),
              width=unit(7,"cm"),
              height=unit(4,"cm"),
              #right_annotation = ra,
              #bottom_annotation = ca,
              #show_row_names = F,
              #show_column_names = F,
              heatmap_legend_param = list(grid_width=unit(3,"mm"),title="# of neighbors (PDAC - PanIN)",title_gp=gpar(fontsize=7),labels_gp=gpar(fontsize=6)),
              name = "# of neighbors")
pdf("../results/Neighbors.pdf",width=7,height=4);draw(hp, merge_legend=T);dev.off()