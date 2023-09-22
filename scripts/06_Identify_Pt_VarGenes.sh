#!/bin/bash
#$ -N PanIN_Pt_VarGenes
#$ -cwd
#$ -pe local 4
#$ -l h_vmem=10G,mem_free=10G,h_fsize=40G
#$ -m e
#$ -M jmitch81@jhmi.edu

module load conda_R/4.1.x
module list

R_SCRIPT="scripts/06_Identify_Pt_VarGenes.R"

R CMD BATCH $R_SCRIPT

exit 0