#!/usr/bin/env Rscript
# 03_cell_QC_and_annotation.R
#
# For one capture: creates a Seurat object from the cellranger output, adds the
# Demuxalot donor assignments, filters low-quality cells and doublets, assigns cell
# types by reference mapping to the PBMC atlas from Hao et al. (2021), and saves the
# QC'd object.
#
# Usage: Rscript 03_cell_QC_and_annotation.R NI-b1
#
# Run with R 4.3.1, Seurat 5.1.0, SeuratObject 5.0.2, SeuratDisk 0.0.0.9021,
# dplyr 1.1.4, ggplot2 3.5.1, patchwork 1.2.0.
# Seurat v5 or later is required since later steps rely on its layer structure.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---------------- User settings: edit these paths ----------------
align_dir      <- "/path/to/alignment"        # OUTPUT_DIR used in 01_alignment.sh
demux_dir      <- "/path/to/demultiplex"      # OUTPUT_DIR used in 02_demultiplexing.sh
reference_file <- "/path/to/multi.h5seurat"   # Hao et al. (2021) PBMC reference
output_dir     <- "/path/to/single-cell_QC"   # QC'd Seurat objects are saved here
plot_dir       <- file.path(output_dir, "plots")

# QC thresholds
min_features   <- 200
max_features   <- 6000
max_percent_mt <- 10
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 03_cell_QC_and_annotation.R <capture>, e.g. NI-b1")
}
capture     <- args[1]                                  # e.g. NI-b1
condition   <- sub("-b[0-9]+$", "", capture)            # e.g. NI
batch       <- as.integer(sub("^.*-b", "", capture))    # e.g. 1
sample_name <- paste0("b", batch, "_", condition)       # e.g. b1_NI, used in file names

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


##########################################
##        CREATE SEURAT OBJECT          ##
##########################################
count_dir <- file.path(align_dir, capture, capture, "outs", "per_sample_outs",
                       capture, "count", "sample_filtered_feature_bc_matrix")
counts <- Read10X(data.dir = count_dir)

# The influenza A gene "NA" (neuraminidase) is read into R as a missing value,
# which would cause it to be dropped. Rename it to "IAV_NA" so it is kept.
na_idx <- which(is.na(rownames(counts)))
if (length(na_idx) > 1) {
  warning("More than one NA gene name found in ", count_dir,
          ": check features.tsv.gz for other missing gene names.")
}
rownames(counts)[na_idx] <- "IAV_NA"

seurat <- CreateSeuratObject(counts = counts, project = paste0(condition, "_b", batch))
rm(counts)


##########################################
##      ADD DEMUXALOT ASSIGNMENTS       ##
##########################################
assignments <- read.delim(gzfile(file.path(demux_dir, capture, "assignments_refined.tsv.gz")))
colnames(assignments) <- c("barcode", "donor")

# Demuxalot labels doublets as "donorA+donorB"; relabel these as "doublet"
assignments$donor <- gsub(".*\\+.*", "doublet", assignments$donor)

# Match assignments to cells by barcode
seurat$donor_id <- assignments$donor[match(colnames(seurat), assignments$barcode)]


##########################################
##                QC                    ##
##########################################
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Violin plots of genes detected, UMI counts, and mitochondrial percentage
plot_qc <- function(obj, stage) {
  plots <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(feature) {
    VlnPlot(obj, features = feature) + theme(legend.position = "none")
  })
  combined <- wrap_plots(plots, ncol = 3) +
    plot_annotation(paste0("batch", batch, " ", condition, " ", stage, "-QC (",
                           ncol(obj), " cells)"))
  ggsave(file.path(plot_dir, paste0(stage, "filter_", sample_name, ".png")),
         combined, width = 10, height = 6, dpi = 600)
}

plot_qc(seurat, "pre")

seurat <- subset(seurat, subset = nFeature_RNA > min_features &
                                  nFeature_RNA < max_features &
                                  percent.mt < max_percent_mt &
                                  donor_id != "doublet")

plot_qc(seurat, "post")


##########################################
##        CELL TYPE ASSIGNMENT          ##
##########################################
# Cells are assigned at three resolutions of the reference annotation:
# celltype.l1 (lowest), celltype.l2, and celltype.l3 (highest).
reference <- LoadH5Seurat(reference_file)

anchors <- FindTransferAnchors(reference = reference,
                               query = seurat,
                               dims = 1:30,
                               normalization.method = "SCT",
                               reference.reduction = "spca")

for (level in names(celltype_levels)) {
  # Transfer labels from the reference and add them to the metadata
  predictions <- TransferData(anchorset = anchors,
                              refdata = reference[[level, drop = TRUE]],
                              dims = 1:30)
  seurat <- AddMetaData(seurat,
                        metadata = setNames(predictions$predicted.id, rownames(predictions)),
                        col.name = level)
}


##########################################
##                SAVE                  ##
##########################################
saveRDS(seurat, file.path(output_dir, paste0(sample_name, "_seurat.rds")))
message("Finished ", capture, ": ", ncol(seurat), " cells retained.")
