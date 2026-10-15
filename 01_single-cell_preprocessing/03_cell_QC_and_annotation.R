#!/usr/bin/env Rscript
# 03_cell_QC_and_annotation.R
#
# For one capture: creates a Seurat object from the Cell Ranger output, adds the
# Demuxalot donor assignments, filters low-quality cells and doublets, assigns cell
# types by reference mapping to the PBMC atlas from Hao, Hao et al., and saves the
# QC'd object along with pre- and post-QC plots.
#
# Usage: Rscript 03_cell_QC_and_annotation.R NI-b1
#
# Run with R 4.3.1, Seurat 5.1.0, SeuratObject 5.0.2, SeuratDisk 0.0.0.9021,
# dplyr 1.1.4, ggplot2 3.5.1, patchwork 1.2.0.
# Seurat v5 or later is required, as later steps rely on its layer structure.

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
reference_file <- "/path/to/multi.h5seurat"   # Hao, Hao et al. PBMC reference
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
##       CREATE SEURAT OBJECT           ##
##########################################
count_dir <- file.path(align_dir, capture, capture, "outs", "per_sample_outs",
                       capture, "count", "sample_filtered_feature_bc_matrix")
counts <- Read10X(data.dir = count_dir)

# Note: the influenza A gene "NA" is read into R as a missing value and dropped
# here; its counts are added back to the pseudobulk data in 06_pbmc_integration_pseudobulk.R.

seurat <- CreateSeuratObject(counts = counts, project = paste0(condition, "_b", batch))
rm(counts)


##########################################
##     ADD DEMUXALOT ASSIGNMENTS        ##
##########################################
assignments <- read.delim(gzfile(file.path(demux_dir, capture, "assignments_refined.tsv.gz")))
colnames(assignments) <- c("barcode", "donor")

# Demuxalot labels doublets as "donorA+donorB"; relabel these as "doublet"
assignments$donor <- gsub(".*\\+.*", "doublet", assignments$donor)

# Match assignments to cells by barcode
seurat$donor_id <- assignments$donor[match(colnames(seurat), assignments$barcode)]


##########################################
##               QC                     ##
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
##      CELL TYPE ASSIGNMENT            ##
##########################################
# Cells are assigned at three resolutions of the reference annotation:
# celltype.l1 (lowest), celltype.l2, and celltype.l3 (highest).
reference <- LoadH5Seurat(reference_file)

anchors <- FindTransferAnchors(reference = reference,
                               query = seurat,
                               dims = 1:30,
                               normalization.method = "SCT",
                               reference.reduction = "spca")

# Order in which cell types are displayed in the distribution plots
celltype_levels <- list(
  celltype.l1 = c("CD4 T", "CD8 T", "NK", "B", "Mono", "other T", "DC", "other"),
  celltype.l2 = c("CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "CD4 CTL",
                  "Treg", "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM", "NK",
                  "NK Proliferating", "NK_CD56bright", "B naive", "B intermediate",
                  "B memory", "Plasmablast", "CD14 Mono", "CD16 Mono", "gdT", "MAIT",
                  "dnT", "cDC1", "cDC2", "ASDC", "pDC", "HSPC", "ILC", "Eryth",
                  "Platelet", "other"),
  celltype.l3 = c("CD4 Naive", "CD4 Proliferating", "CD4 TCM_1", "CD4 TCM_2", "CD4 TCM_3",
                  "CD4 TEM_1", "CD4 TEM_2", "CD4 TEM_3", "CD4 TEM_4", "CD4 CTL",
                  "Treg Naive", "Treg Memory", "CD8 Naive", "CD8 Naive_2",
                  "CD8 Proliferating", "CD8 TCM_1", "CD8 TCM_2", "CD8 TCM_3", "CD8 TEM_1",
                  "CD8 TEM_2", "CD8 TEM_3", "CD8 TEM_4", "CD8 TEM_5", "CD8 TEM_6", "NK_1",
                  "NK_2", "NK_3", "NK_4", "NK Proliferating", "NK_CD56bright",
                  "B naive kappa", "B naive lambda", "B intermediate kappa",
                  "B intermediate lambda", "B memory kappa", "B memory lambda",
                  "Plasmablast", "CD14 Mono", "CD16 Mono", "gdT_1", "gdT_2", "gdT_3",
                  "gdT_4", "MAIT", "dnT_1", "dnT_2", "cDC1", "cDC2_1", "cDC2_2",
                  "ASDC_mDC", "ASDC_pDC", "pDC", "HSPC", "ILC", "Eryth", "Platelet",
                  "Plasma")
)
plot_widths <- c(celltype.l1 = 10, celltype.l2 = 10, celltype.l3 = 15)

for (level in names(celltype_levels)) {
  # Transfer labels from the reference and add them to the metadata
  predictions <- TransferData(anchorset = anchors,
                              refdata = reference[[level, drop = TRUE]],
                              dims = 1:30)
  seurat <- AddMetaData(seurat,
                        metadata = setNames(predictions$predicted.id, rownames(predictions)),
                        col.name = level)

  # Plot the proportion of cells assigned to each cell type
  plot_df <- data.frame(celltype = seurat@meta.data[[level]]) %>%
    count(celltype, name = "count") %>%
    mutate(percent = count / sum(count),
           celltype = factor(celltype, levels = celltype_levels[[level]]))

  p <- ggplot(plot_df, aes(x = celltype, y = percent)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggtitle(paste0("batch", batch, " ", condition, " ", level)) +
    labs(x = "cell type", y = "proportion") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12))

  ggsave(file.path(plot_dir, paste0(gsub(".", "", level, fixed = TRUE),
                                    "_distribution_", sample_name, ".png")),
         p, width = plot_widths[[level]], height = 6, dpi = 600)
}


##########################################
##              SAVE                    ##
##########################################
saveRDS(seurat, file.path(output_dir, paste0(sample_name, "_seurat.rds")))
message("Finished ", capture, ": ", ncol(seurat), " cells retained.")
