#!/usr/bin/env Rscript
# 07_capture_corrected_expression_matrices.R
# For each of the five main cell types: keeps protein-coding genes and samples in
# the metadata, keeps genes highly expressed in at least one cell type, and
# regresses out capture effects while keeping condition effects. Saves corrected
# expression and voom weights for all samples and for two filtered sample sets.
# Usage: Rscript 07_capture_corrected_expression_matrices.R
# Run with R 4.3.3, limma 3.58.1, edgeR 4.0.16.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})

# ---------------- User settings: edit these paths ----------------
# Paths are relative to this folder; see data/README.md
counts_file     <- "../data/integration/raw_counts_list_unfiltered.rds"
metadata_dir    <- "../data/metadata"
output_dir      <- "../data/counts_and_weights"
ensembl_release <- 115                                       # Ensembl release used for gene biotypes
# Minimum median log-CPM for a gene to count as highly expressed in each cell type
median_thresholds <- c(NK = 2.5, CD8T = 2.5, monocyte = 1.5, B = 2, CD4T = 2)
# ------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
counts_list    <- readRDS(counts_file)
meta_all       <- readRDS(file.path(metadata_dir, "meta_list_unfiltered.rds"))
meta_missingNI  <- readRDS(file.path(metadata_dir, "meta_list_filtered_missingNI.rds"))
meta_missingBMI <- readRDS(file.path(metadata_dir, "meta_list_filtered_missingBMI.rds"))
celltypes <- setNames(names(meta_all), names(meta_all))

# Normalize with voom, modeling condition and capture (batch within condition)
run_voom <- function(counts, meta) {
  dge <- calcNormFactors(DGEList(counts = counts))
  design <- model.matrix(~ 0 + condition + condition:batch, data = meta)
  design <- design[, colSums(design != 0) > 0]
  list(v = voom(dge, design, plot = FALSE), design = design)
}


# Get human protein-coding genes from the Ensembl gene annotation (GTF)
gtf_file <- file.path(tempdir(), "ensembl.gtf.gz")
options(timeout = 1200)   # allow up to 20 minutes for the download
download.file(sprintf("https://ftp.ensembl.org/pub/release-%d/gtf/homo_sapiens/Homo_sapiens.GRCh38.%d.gtf.gz",
                      ensembl_release, ensembl_release), gtf_file, mode = "wb")
gtf <- readLines(gzfile(gtf_file))
gtf <- gtf[grepl("\tgene\t", gtf) & grepl('gene_biotype "protein_coding"', gtf, fixed = TRUE) &
           grepl("gene_name", gtf, fixed = TRUE)]
protein_coding <- unique(sub('.*gene_name "([^"]+)".*', "\\1", gtf))

# Keep protein-coding genes and samples in the metadata (samples with at least 10 cells making up pseudobulk)
counts_list <- lapply(celltypes, function(ct) {
  counts <- as.matrix(counts_list[[ct]])
  counts <- counts[rownames(counts) %in% protein_coding, , drop = FALSE]
  counts <- counts[, colnames(counts) %in% rownames(meta_all[[ct]]), drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(meta_all[[ct]])))   # samples must match metadata order
  counts
})

# Keep genes highly expressed in at least one cell type (union across cell types)
genes_keep <- unique(unlist(lapply(celltypes, function(ct) {
  v <- run_voom(counts_list[[ct]], meta_all[[ct]])$v
  rownames(v$E)[apply(v$E, 1, median) > median_thresholds[[ct]]]
})))

# Regress out capture effects, keeping condition effects
corrected <- list()
weights <- list()
for (ct in celltypes) {
  res <- run_voom(counts_list[[ct]][genes_keep, ], meta_all[[ct]])
  fit <- eBayes(lmFit(res$v, res$design))
  condition_cols <- grep("^condition(NI|LPS|IAV|Tcruzi)$", colnames(res$design), value = TRUE)
  condition_effect <- fit$coefficients[, condition_cols] %*% t(res$design[, condition_cols])
  corrected[[ct]] <- residuals.MArrayLM(fit, res$v) + condition_effect
  weights[[ct]] <- res$v$weights
  dimnames(weights[[ct]]) <- dimnames(corrected[[ct]])
}


# Save for all samples, and for samples from donors with NI data or with BMI data
save_set <- function(suffix, meta_list) {
  keep <- function(m, ct) m[, colnames(m) %in% rownames(meta_list[[ct]]), drop = FALSE]
  saveRDS(lapply(celltypes, function(ct) keep(corrected[[ct]], ct)),
          file.path(output_dir, paste0("corrected_counts_list_", suffix, ".rds")))
  saveRDS(lapply(celltypes, function(ct) keep(weights[[ct]], ct)),
          file.path(output_dir, paste0("weights_list_", suffix, ".rds")))
}
save_set("unfiltered", meta_all)
save_set("filtered_missingNI", meta_missingNI)
save_set("filtered_missingBMI", meta_missingBMI)
message("Finished")
