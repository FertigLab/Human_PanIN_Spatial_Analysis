#!/bin/bash
#SBATCH --job-name=TLS_visium_aggregation
#SBATCH --time=24:00:00
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128G
#SBATCH --mail-type=end
#SBATCH --mail-user=jmitch81@jhmi.edu
#SBATCH -o ./slurm-%A_%a.out

module load seurat/4.1.1

module list

R_SCRIPT="scripts/16_A_TLS_Selection_for_CoGAPS.R"

R CMD BATCH $R_SCRIPT

echo "Finished with job $SLURM_JOBID"
