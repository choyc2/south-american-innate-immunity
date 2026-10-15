#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --time=24:00:00
#SBATCH --mem=250G
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Runs 04_merge.R to combine all captures into a single Seurat object.
# Edit the paths at the top of the R script before running.
#
# Usage (from this folder, after 03_cell_QC_and_annotation.sh has finished
# for all captures):
#   mkdir -p logs
#   sbatch 04_merge.sh

# Modules used on our cluster; adjust these to match your system
module load gsl/2.7
module load gcc/10.2.0
module load hdf5/1.12.0
module load libxc/5.1.7+intelmpi-2019.up7+intel-19.1.1
module load R/4.3.1

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"   # so 04_merge.R finds captures.txt in this folder
Rscript 04_merge.R
