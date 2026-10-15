#!/bin/bash
#SBATCH --job-name=cell_QC
#SBATCH --time=8:00:00
#SBATCH --mem=200G
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Runs 03_cell_QC_and_annotation.R for one capture.
# Edit the paths at the top of the R script before running.
#
# Usage (from this folder, after 02_demultiplexing.sh has finished for the capture):
#   mkdir -p logs
#   sbatch --job-name=qc_NI-b1 03_cell_QC_and_annotation.sh NI-b1

# Modules used on our cluster; adjust these to match your system
module load gsl/2.7
module load gcc/10.2.0
module load hdf5/1.12.0
module load libxc/5.1.7+intelmpi-2019.up7+intel-19.1.1
module load R/4.3.1

set -euo pipefail

CAPTURE=${1:?"Usage: sbatch 03_cell_QC_and_annotation.sh <capture>, e.g. NI-b1"}

Rscript "${SLURM_SUBMIT_DIR:-.}/03_cell_QC_and_annotation.R" "${CAPTURE}"
