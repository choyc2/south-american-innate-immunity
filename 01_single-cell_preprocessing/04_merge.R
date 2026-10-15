#!/usr/bin/env Rscript
# 04_merge.R
#
# Merges the QC'd Seurat objects from all captures into a single, unnormalized
# Seurat object. Captures are first merged within each condition (one layer per
# batch, then joined), and the four conditions are then merged together.
# Condition, batch, and 10x chemistry are added to the metadata, and donors
# not included in this study are removed.
#
# Usage (from this folder, after 03_cell_QC_and_annotation.R has finished
# for all captures):
#   Rscript 04_merge.R
#
# Run with R 4.3.1, Seurat 5.1.0, SeuratObject 5.0.2, dplyr 1.1.4, stringr 1.5.1.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

# ---------------- User settings: edit these paths ----------------
qc_dir        <- "/path/to/single-cell_QC"   # output_dir used in 03_cell_QC_and_annotation.R
captures_file <- "captures.txt"
output_file   <- "/path/to/merged_seurat_unnormalized.rds"
# ------------------------------------------------------------------

captures <- readLines(captures_file)
captures <- captures[captures != ""]
capture_info <- data.frame(capture   = captures,
                           condition = sub("-b[0-9]+$", "", captures),
                           batch     = as.integer(sub("^.*-b", "", captures)))
capture_info <- capture_info[order(capture_info$batch), ]


##########################################
##     MERGE WITHIN EACH CONDITION      ##
##########################################
conditions <- c("NI", "LPS", "IAV", "Tcruzi")
condition_objects <- list()

for (cond in conditions) {
  batches <- capture_info$batch[capture_info$condition == cond]
  sample_names <- paste0("b", batches, "_", cond)   # e.g. b1_NI; also used as cell name prefixes
  message(cond, ": merging ", length(sample_names), " batches")

  objects <- lapply(sample_names, function(s) {
    readRDS(file.path(qc_dir, paste0(s, "_seurat.rds")))
  })

  merged <- merge(x = objects[[1]],
                  y = objects[-1],
                  add.cell.ids = sample_names,
                  merge.data = FALSE,   # only merge counts
                  merge.dr = FALSE,     # don't merge dimensional reductions
                  project = cond)

  condition_objects[[cond]] <- JoinLayers(merged)
  rm(objects, merged)
  gc()
}


##########################################
##        MERGE ALL CONDITIONS          ##
##########################################
merged_seurat_unnormalized <- merge(x = condition_objects[["NI"]],
                                    y = condition_objects[c("LPS", "IAV", "Tcruzi")],
                                    merge.data = FALSE,
                                    merge.dr = FALSE,
                                    project = "CZI_innate")
rm(condition_objects)
gc()


##########################################
##          EDIT THE METADATA           ##
##########################################
meta <- merged_seurat_unnormalized@meta.data

# Condition, from orig.ident (e.g. NI_b1 -> NI)
merged_seurat_unnormalized$condition <- case_when(
  grepl("NI", meta$orig.ident)     ~ "NI",
  grepl("LPS", meta$orig.ident)    ~ "LPS",
  grepl("IAV", meta$orig.ident)    ~ "IAV",
  grepl("Tcruzi", meta$orig.ident) ~ "Tcruzi",
  TRUE                             ~ NA_character_
)

# Batch number, from orig.ident (e.g. NI_b1 -> 1)
merged_seurat_unnormalized$batch <- as.numeric(str_extract(meta$orig.ident, "\\d+"))

# 10x chemistry: batches 1-9 used Next GEM, batches 10-16 used GEM-X
merged_seurat_unnormalized$technology <- case_when(
  merged_seurat_unnormalized$batch >= 1  & merged_seurat_unnormalized$batch <= 9  ~ "NEXT_GEM",
  merged_seurat_unnormalized$batch >= 10 & merged_seurat_unnormalized$batch <= 16 ~ "GEM-X",
  TRUE                                                                          ~ NA_character_
)

# Reorder metadata columns
merged_seurat_unnormalized@meta.data <- merged_seurat_unnormalized@meta.data[
  c("orig.ident", "condition", "batch", "donor_id", "nCount_RNA", "nFeature_RNA",
    "percent.mt", "technology", "celltype.l1", "celltype.l2", "celltype.l3")
]
rm(meta)


##########################################
##       REMOVE EXCLUDED DONORS         ##
##########################################
# RIH donors (individuals of European ancestry living in the US) are not
# included in this study
cells_keep <- colnames(merged_seurat_unnormalized)[
  !grepl("^RIH", merged_seurat_unnormalized$donor_id)
]
merged_seurat_unnormalized <- subset(merged_seurat_unnormalized, cells = cells_keep)

# Identify donors sequenced in more than one batch (each donor is expected to
# appear in at most four captures, one per condition)
duplicate_donors <- merged_seurat_unnormalized@meta.data %>%
  distinct(donor_id, orig.ident) %>%
  count(donor_id) %>%
  filter(n > 4)
print(duplicate_donors)
# CHI126 and MAP059 were each sequenced in batches 7 and 10. The batch 10
# samples were kept because they used the newer GEM-X chemistry.
cells_keep <- colnames(merged_seurat_unnormalized)[
  !(merged_seurat_unnormalized$batch == 7 &
    merged_seurat_unnormalized$donor_id %in% c("CHI126", "MAP059"))
]
merged_seurat_unnormalized <- subset(merged_seurat_unnormalized, cells = cells_keep)


##########################################
##                SAVE                  ##
##########################################
saveRDS(merged_seurat_unnormalized, output_file)
message("Finished: ", ncol(merged_seurat_unnormalized), " cells from ",
        length(unique(merged_seurat_unnormalized$donor_id)), " donors.")
