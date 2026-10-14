#!/bin/bash
#SBATCH --job-name=demuxalot
#SBATCH --time=10:00:00
#SBATCH --mem=80G
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Assigns cells in one capture to individual donors using Demuxalot
# (run through the Demuxafy Singularity container), based on donor genotypes.
#
# Usage (from this folder, after 01_alignment.sh has finished for the capture):
#   mkdir -p logs
#   sbatch --job-name=demux_NI-b1 02_demultiplexing.sh NI-b1
#
# Capture names must follow the format <condition>-b<batch number>, e.g. NI-b1.
# The batch number is used to select the donor list and genotype VCF for that batch.

module load singularity

set -euo pipefail

CAPTURE=${1:?"Usage: sbatch 02_demultiplexing.sh <capture>, e.g. NI-b1"}
BATCH=${CAPTURE##*-b}    # e.g. NI-b1 -> 1

# ---------------- User settings: edit these paths ----------------
DEMUXAFY_SIF=/path/to/Demuxafy.sif     # Demuxafy Singularity image
BIND_PATH=/path/to/project              # parent folder containing all input and output paths below
ALIGN_DIR=/path/to/alignment            # OUTPUT_DIR used in 01_alignment.sh
GENOTYPE_DIR=/path/to/genotypes         # per-batch donor lists and VCFs (see README)
OUTPUT_DIR=/path/to/demultiplex         # one subfolder per capture is created here
# ------------------------------------------------------------------

COUNT_DIR="${ALIGN_DIR}/${CAPTURE}/${CAPTURE}/outs/per_sample_outs/${CAPTURE}/count"
mkdir -p "${OUTPUT_DIR}/${CAPTURE}"

singularity exec --bind "${BIND_PATH}" "${DEMUXAFY_SIF}" Demuxalot.py \
    -a "${COUNT_DIR}/sample_alignments.bam" \
    -b "${COUNT_DIR}/sample_filtered_feature_bc_matrix/barcodes.tsv.gz" \
    -n "${GENOTYPE_DIR}/batch${BATCH}_inds.txt" \
    -v "${GENOTYPE_DIR}/batch${BATCH}.MAF0.05.GRCh38_v2.vcf.gz" \
    -o "${OUTPUT_DIR}/${CAPTURE}" \
    -r True    # refine donor genotypes using the single-cell data
