#!/usr/bin/env Rscript
# 06_integration_pseudobulk.R
# Integrates all cells kept by 05_celltype_annotation.R, then creates pseudobulk
# counts per donor and capture for the five main cell types (CD4T, CD8T, NK,
# monocyte, B).
# Usage: Rscript 06_integration_pseudobulk.R
# Run with R 4.3.1, Seurat 5.1.0, dplyr 1.1.4, stringr 1.5.1.

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(dplyr)
  library(stringr)
})

# ---------------- User settings: edit these paths ----------------
merged_file     <- "/path/to/merged_seurat_unnormalized.rds"   # output of 04_merge.R
integration_dir <- "../data/integration"                       # relative to this folder; see data/README.md
align_dir       <- "/path/to/alignment"                        # OUTPUT_DIR in 01_alignment.sh
captures_file   <- "captures.txt"
# ------------------------------------------------------------------
dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)


##########################################
##          FINAL INTEGRATION           ##
##########################################
celltypes_all <- c("B", "CD4T", "CD8T", "DC", "monocyte", "NK", "otherT", "other")
keep <- do.call(rbind, lapply(celltypes_all, function(ct) {
  read.csv(file.path(integration_dir, paste0(ct, "_keep.csv")))
}))

seurat <- readRDS(merged_file)
seurat <- subset(seurat, cells = intersect(Cells(seurat), keep$cell_barcode))
idx <- match(colnames(seurat), keep$cell_barcode)
seurat$celltype.l1_manual <- keep$celltype.l1_manual[idx]
seurat$celltype.l2_manual <- keep$celltype.l2_manual[idx]

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
seurat <- RunUMAP(seurat, reduction = "integrated.rpca", dims = 1:60,
                  reduction.name = "rna.umap", reduction.key = "UMAP_", min.dist = 0.3)

seurat@meta.data <- seurat@meta.data[c("orig.ident", "condition", "batch", "donor_id",
                                       "nCount_RNA", "nFeature_RNA", "percent.mt",
                                       "technology", "celltype.l1", "celltype.l2",
                                       "celltype.l3", "celltype.l1_manual",
                                       "celltype.l2_manual", "integrated_clusters")]

metadata <- cbind(seurat@meta.data, Embeddings(seurat, reduction = "rna.umap"))
saveRDS(metadata, file.path(integration_dir, "final_pbmc_seurat_south_america_innate_metadata.rds"))


##########################################
##              PSEUDOBULK              ##
##########################################
# Number of cells per donor, capture, and cell type. AggregateExpression
# replaces "_" with "-" in sample names, so match that here.
cell_counts <- metadata %>%
  count(orig.ident, donor_id, celltype.l1_manual, name = "count") %>%
  mutate(orig.ident = str_replace_all(orig.ident, "_", "-"))

# Sum counts per donor and capture, and add batch, condition, and cell numbers
make_pseudobulk <- function(obj, counts_df, label) {
  pb <- AggregateExpression(obj, group.by = c("orig.ident", "donor_id"), return.seurat = TRUE)
  md <- pb@meta.data
  md$batch <- as.numeric(str_extract(md$orig.ident, "\\d+"))
  md$condition <- str_extract(md$orig.ident, "NI|LPS|IAV|Tcruzi")
  md$orig.ident <- str_remove(md$orig.ident, "_.*$")
  md <- left_join(md, counts_df, by = c("orig.ident", "donor_id"))
  md$celltype <- label
  rownames(md) <- paste(md$orig.ident, md$donor_id, md$celltype, sep = "-")
  pb@meta.data <- md
  pb
}

main_celltypes <- c("CD4T", "CD8T", "NK", "monocyte", "B")
pseudo_list <- list()
for (ct in main_celltypes) {
  counts_ct <- cell_counts %>% filter(celltype.l1_manual == ct) %>%
    select(orig.ident, donor_id, count)
  pseudo_list[[ct]] <- make_pseudobulk(subset(seurat, subset = celltype.l1_manual == ct),
                                       counts_ct, ct)
}


##########################################
##     ADD THE INFLUENZA "NA" GENE      ##
##########################################
# The influenza A gene "NA" is read into R as a missing value and dropped when the
# Seurat objects are created. Its counts are read directly from the Cell Ranger
# output and added to the cell type pseudobulk count matrices as "IAV_NA".
captures <- readLines(captures_file)
na_counts <- unlist(unname(lapply(captures[captures != ""], function(capture) {
  mat <- Read10X(file.path(align_dir, capture, capture, "outs", "per_sample_outs",
                           capture, "count", "sample_filtered_feature_bc_matrix"))
  na_idx <- which(is.na(rownames(mat)))
  if (length(na_idx) == 0) return(NULL)
  condition <- sub("-b[0-9]+$", "", capture)
  batch <- sub("^.*-b", "", capture)
  setNames(mat[na_idx, ], paste0("b", batch, "_", condition, "_", colnames(mat)))
})))

na_pseudobulk <- metadata %>%
  mutate(na_count = coalesce(na_counts[rownames(metadata)], 0),
         orig.ident = str_replace_all(orig.ident, "_", "-")) %>%
  filter(celltype.l1_manual %in% main_celltypes) %>%
  group_by(orig.ident, donor_id, celltype.l1_manual) %>%
  summarise(na_count = sum(na_count), .groups = "drop") %>%
  mutate(sample = paste(orig.ident, donor_id, celltype.l1_manual, sep = "-"))

# Count matrices per cell type, with samples named <capture>-<donor>-<cell type>
raw_counts_list_unfiltered <- lapply(pseudo_list, function(pb) {
  m <- LayerData(pb, layer = "counts")
  colnames(m) <- rownames(pb@meta.data)
  m
})
for (ct in main_celltypes) {
  m <- raw_counts_list_unfiltered[[ct]]
  na_vec <- setNames(rep(0, ncol(m)), colnames(m))
  sub <- na_pseudobulk[na_pseudobulk$celltype.l1_manual == ct, ]
  matched <- intersect(colnames(m), sub$sample)
  na_vec[matched] <- sub$na_count[match(matched, sub$sample)]
  raw_counts_list_unfiltered[[ct]] <- rbind(m, IAV_NA = na_vec)
}
saveRDS(raw_counts_list_unfiltered, file.path(integration_dir, "raw_counts_list_unfiltered.rds"))
message("Finished")
