#!/bin/bash
#SBATCH --job-name=18_B_TLS_SCENIC
#SBATCH --time=72:00:00
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=150G
#SBATCH --mail-type=end
#SBATCH --mail-user=jmitch81@jhmi.edu
#SBATCH -o ./reports/slurm-%A_%a.out

module load seurat/4.1.1
module list

DATA_DIR="processed_data/18_A_TLS_domino_preprocessing"
RESULT_DIR="processed_data/18_B_TLS_SCENIC"
ID="TLS_neighbors"

# mkdir ${RESULT_DIR}

# input data
LOOM_FILE="${DATA_DIR}/TLS_neighbors_counts.loom"
TF_REFERENCE="reference/allTFs_hg38.txt"
TSS_REFERENCE_1="reference/HG38/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
TSS_REFERENCE_2="reference/HG38/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
MOTIFS_REFERENCE="reference/HG38/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

# generated results
GRN="${RESULT_DIR}/adj_${ID}.tsv"
REGULON="${RESULT_DIR}/regulons_${ID}.csv"
AUC="${RESULT_DIR}/auc_${ID}.csv"

# calculate gene regulatory network
singularity exec ~/venv/domino/aertslab-pyscenic-0.11.0.sif pyscenic grn \
    ${LOOM_FILE} \
    ${TF_REFERENCE} \
    -o $GRN \
    --num_workers 8

echo "${ID} Step 1: gene regulatory network calculation complete"

# calculate regulons
singularity exec ~/venv/domino/aertslab-pyscenic-0.11.0.sif pyscenic ctx \
    ${GRN} \
    ${TSS_REFERENCE_1} \
    ${TSS_REFERENCE_2} \
    --annotations_fname ${MOTIFS_REFERENCE} \
    --expression_mtx_fname "${LOOM_FILE}" \
    --mode "dask_multiprocessing" \
    --output "${REGULON}" \
    --num_workers 8
    
echo "${ID} Step 2: regulons calculation complete"

# AUC / TF activity calculation

singularity exec ~/venv/domino/aertslab-pyscenic-0.11.0.sif pyscenic aucell \
  	"${LOOM_FILE}" \
  	"${REGULON}" \
  	-o "${AUC}"
    
echo "${PT_ID} Step 3: AUC matrix calculation complete"

exit 0