#!/bin/bash
#SBATCH --job-name=06_integration_pseudobulk
#SBATCH --time=30:00:00
#SBATCH --mem=400G
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Usage (from this folder): sbatch 06_pbmc_integration_pseudobulk.sh

module load gsl/2.7
module load gcc/10.2.0
module load hdf5/1.12.0
module load libxc/5.1.7+intelmpi-2019.up7+intel-19.1.1
module load R/4.3.1

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
Rscript 06_pbmc_integration_pseudobulk.R
