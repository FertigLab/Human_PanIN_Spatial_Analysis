# 16_B_12_Pattern_TLS_CoGAPS

library(CoGAPS)
library(Matrix)

n_patterns <- 12

# create results directory
result_dir <- paste0("processed_data/16_B_", n_patterns, "_Pattern_TLS_CoGAPS")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

# load log-transformed counts matrix
mat <- read.csv(paste0("processed_data/",
                       "16_A_TLS_Selection_for_CoGAPS/",
                       "TLS_visium_logCounts.csv"),
                header = TRUE)[,-1]
dim(mat)

# load parameters
params <- readRDS(paste0("processed_data/",
                         "16_A_TLS_Selection_for_CoGAPS/",
                         "CoGAPS_parameters_", n_patterns, "_patterns.rds"))
print(params)
cogaps_result <- CoGAPS(mat, params = params)

saveRDS(cogaps_result,
        file = paste0(result_dir, "/cogaps_result_", n_patterns, "_patterns.rds"))

rm(list = ls())