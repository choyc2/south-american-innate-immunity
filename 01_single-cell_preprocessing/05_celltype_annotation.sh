#!/bin/bash
#SBATCH --job-name=05_celltype_annotation
#SBATCH --time=20:00:00
#SBATCH --mem=300G
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Usage (from this folder): sbatch --job-name=annot_monocyte 05_celltype_annotation.sh monocyte

module load gsl/2.7
module load gcc/10.2.0
module load hdf5/1.12.0
module load libxc/5.1.7+intelmpi-2019.up7+intel-19.1.1
module load R/4.3.1

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
Rscript 05_celltype_annotation.R "${1:?Give a cell type, e.g. monocyte}"
