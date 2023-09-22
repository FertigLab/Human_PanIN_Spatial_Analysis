#!/bin/bash
#$ -N PanIN_Vis_Var_gene
#$ -cwd
#$ -pe local 4
#$ -l h_vmem=12G,mem_free=12G,h_fsize=40G
#$ -m e
#$ -M jmitch81@jhmi.edu

module load conda_R/4.1.x
module list

R_SCRIPT="scripts/01_Read_Segments_Normalize_and_Scale.R"

R CMD BATCH $R_SCRIPT

exit 0
