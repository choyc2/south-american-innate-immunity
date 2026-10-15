#!/usr/bin/env Rscript
# 05_celltype_annotation.R
# For one broad cell type: integrates across conditions, clusters, and applies the
# manual annotation decisions below. Writes the cells kept and their final labels.
# Usage: Rscript 05_celltype_annotation.R monocyte
# Run with R 4.3.1, Seurat 5.1.0.

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
})

# ---------------- User settings: edit these paths ----------------
merged_file     <- "/path/to/merged_seurat_unnormalized.rds"   # output of 04_merge.R
integration_dir <- "/path/to/integration"
# ------------------------------------------------------------------

# Manual annotation decisions are based on inspection of marker genes (see Methods
# section). A cell is kept if its celltype.l2 and celltype.l3 assignments are valid
# subsets of the cell's celltype.l1 assignment.
tcm_tem <- function(type, n_tcm, n_tem) {
  c(setNames(rep(paste0(type, "_TCM"), n_tcm), paste0(type, " TCM_", 1:n_tcm)),
    setNames(rep(paste0(type, "_TEM"), n_tem), paste0(type, " TEM_", 1:n_tem)))
}
decisions <- list(
  monocyte = list(l1 = "Mono", label = "monocyte",
                  l3_to_label = c("CD14 Mono" = "CD14_monocyte", "CD16 Mono" = "CD16_monocyte"),
                  l2 = c("CD14 Mono", "CD16 Mono"),
                  keep_clusters = c(0:3, 5:8, 10:16, 19:21)),
  CD4T = list(l1 = "CD4 T", label = "CD4T",
              l3_to_label = c("CD4 CTL" = "CD4 CTL", "CD4 Naive" = "CD4 Naive",
                              "CD4 Proliferating" = "CD4 Proliferating",
                              tcm_tem("CD4", 3, 4),
                              "Treg Memory" = "Treg", "Treg Naive" = "Treg"),
              l2 = c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "Treg"),
              keep_clusters = c(0:8, 12, 13, 15, 16, 18, 19, 22)),
  CD8T = list(l1 = "CD8 T", label = "CD8T",
              l3_to_label = c("CD8 Naive" = "CD8_naive", "CD8 Naive_2" = "CD8_naive",
                              "CD8 Proliferating" = "CD8_proliferating",
                              tcm_tem("CD8", 3, 6)),
              l2 = c("CD8 Naive", "CD8 TEM", "CD8 TCM", "CD8 Proliferating"),
              keep_clusters = c(0:11, 15:17, 19, 23:25, 27, 29)),
  B = list(l1 = "B", label = "B",
           l3_to_label = c("B naive kappa" = "B_naive", "B naive lambda" = "B_naive",
                           "B intermediate kappa" = "B_intermediate",
                           "B intermediate lambda" = "B_intermediate",
                           "B memory kappa" = "B_memory", "B memory lambda" = "B_memory",
                           "Plasmablast" = "plasmablast", "Plasma" = "plasmablast"),
           l2 = c("B naive", "B intermediate", "B memory", "Plasmablast"),
           keep_clusters = c(0:9, 12, 14:17, 21, 22)),
  NK = list(l1 = "NK", label = "NK",
            l3_to_label = c("NK_1" = "NK_CD56dim", "NK_2" = "NK_CD56dim",
                            "NK_3" = "NK_CD56dim", "NK_4" = "NK_CD56dim",
                            "NK_CD56bright" = "NK_CD56bright",
                            "NK Proliferating" = "NK_proliferating"),
            l2 = c("NK", "NK_CD56bright", "NK Proliferating"),
            keep_clusters = c(0:5, 7, 8, 11, 14, 15, 20)),
  DC = list(l1 = "DC", label = "DC",
            l3_to_label = c("cDC1" = "DC_conventional", "cDC2_1" = "DC_conventional",
                            "cDC2_2" = "DC_conventional", "pDC" = "DC_plasmacytoid",
                            "ASDC_mDC" = "DC_AXL", "ASDC_pDC" = "DC_AXL"),
            l2 = c("cDC1", "cDC2", "pDC", "ASDC"),
            keep_clusters = c(1:6, 9, 11, 13)),
  otherT = list(l1 = "other T", label = "otherT",
                l3_to_label = c("gdT_1" = "gdT", "gdT_2" = "gdT", "gdT_3" = "gdT", "gdT_4" = "gdT",
                                "dnT_1" = "dnT", "dnT_2" = "dnT", "MAIT" = "MAIT"),
                l2 = c("gdT", "MAIT", "dnT"),
                keep_clusters = c(0, 1, 3, 5:7)),
  other = list(l1 = "other", label = "HSPC",
               l3_to_label = c("HSPC" = "HSPC"),
               l2 = "HSPC",
               keep_clusters = 11)
)

celltype <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(celltype) || !celltype %in% names(decisions)) {
  stop("Usage: Rscript 05_celltype_annotation.R <one of: ",
       paste(names(decisions), collapse = ", "), ">")
}
d <- decisions[[celltype]]
target_l1 <- d$l1

# Subset to this cell type
seurat <- readRDS(merged_file)
seurat <- subset(seurat, subset = celltype.l1 == target_l1)

# Normalize, integrate across conditions, and cluster
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat, vars.to.regress = c("batch", "percent.mt", "technology"))
seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = 60)
plan(sequential)
options(future.globals.maxSize = 10 * 1024^3)
seurat <- IntegrateLayers(seurat, method = RPCAIntegration, orig.reduction = "pca",
                          new.reduction = "integrated.rpca",
                          normalization.method = "LogNormalize")
seurat <- FindNeighbors(seurat, dims = 1:30, reduction = "integrated.rpca")
seurat <- FindClusters(seurat, resolution = 0.8, cluster.name = "integrated_clusters")

# Apply the decisions and save the cells kept
meta <- seurat@meta.data
keep <- meta$celltype.l3 %in% names(d$l3_to_label) &
        meta$celltype.l2 %in% d$l2 &
        meta$integrated_clusters %in% as.character(d$keep_clusters)

write.csv(data.frame(cell_barcode       = rownames(meta)[keep],
                     celltype.l1_manual = d$label,
                     celltype.l2_manual = unname(d$l3_to_label[meta$celltype.l3[keep]])),
          file.path(integration_dir, paste0(celltype, "_keep.csv")), row.names = FALSE)
message(celltype, ": kept ", sum(keep), " of ", nrow(meta), " cells")
