#!/bin/bash
#SBATCH --job-name=cellranger_align
#SBATCH --time=30:00:00
#SBATCH --mem=1G
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Aligns reads for one capture with cellranger multi (v9.0.1).
# Low memory is requested because --jobmode=slurm submits the heavy
# alignment steps as separate SLURM jobs.
#
# Usage (from this folder):
#   mkdir -p logs
#   sbatch --job-name=align_IAV-b1 01_alignment.sh IAV-b1

set -euo pipefail

CAPTURE=${1:?"Usage: sbatch 01_alignment.sh <capture>, e.g. IAV-b1"}

# ---------------- User settings: edit these paths ----------------
CELLRANGER=/path/to/cellranger-9.0.1/cellranger
REFERENCE=/path/to/hg38_iav_tcNonEsm      # combined human, IAV, and T. cruzi reference
FASTQ_DIR=/path/to/fastqs                 # contains one subfolder per capture
OUTPUT_DIR=/path/to/alignment             # one subfolder per capture is created here
# ------------------------------------------------------------------

CONFIG_TEMPLATE="${SLURM_SUBMIT_DIR:-.}/configs/${CAPTURE}_config.csv"
mkdir -p "${OUTPUT_DIR}/${CAPTURE}"
CONFIG="${OUTPUT_DIR}/${CAPTURE}/${CAPTURE}_config.csv"

# Fill in the placeholder paths in the capture's config template
sed -e "s|REFERENCE_PATH|${REFERENCE}|g" \
    -e "s|FASTQ_DIR|${FASTQ_DIR}|g" \
    "${CONFIG_TEMPLATE}" > "${CONFIG}"

cd "${OUTPUT_DIR}/${CAPTURE}"

"${CELLRANGER}" multi \
    --id="${CAPTURE}" \
    --csv="${CONFIG}" \
    --jobmode=slurm \
    --maxjobs=32
