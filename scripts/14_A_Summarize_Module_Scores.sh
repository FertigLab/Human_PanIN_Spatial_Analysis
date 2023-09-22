#!/bin/bash
#SBATCH --job-name=mod_score_summary
#SBATCH --time=24:00:00
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --mail-type=end
#SBATCH --mail-user=jmitch81@jhmi.edu
#SBATCH -o ./slurm-%A_%a.out

module load seurat/4.1.1

module list

R_SCRIPT="scripts/14_A_Summarize_Module_Scores.R"

R CMD BATCH $R_SCRIPT

echo "Finished with job $SLURM_JOBID"
