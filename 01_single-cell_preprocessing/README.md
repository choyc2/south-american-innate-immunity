# 01. Single-cell preprocessing

Scripts to process raw single-cell RNA-seq reads into a batch-corrected pseudobulk gene expression matrix. Run the scripts in numbered order.

| Script | Description | Input | Output |
|--------|-------------|-------|--------|
| `01_alignment.sh` | Aligns reads and generates count matrices using cellranger multi (10x Genomics, v9.0.1) | Raw FASTQ files, config files | Gene-by-cell count matrices |
| `02_demultiplexing.sh` | Assigns cells to individuals using genotype data with demuxalot (v0.4.1) | Count matrices, genotype VCF | Cell-to-individual assignments |
| `03_cell_QC.R` | Removes low-quality cells and doublets | Count matrices, assignments | Filtered cells |
| `04_annotation.R` | Annotates cell types | Filtered cells, reference dataset | Annotated cells |
| `05_pseudobulk.R` | Sums counts per individual, cell type, and condition | Annotated cells | Pseudobulk counts |
| `06_batch_correction.R` | Normalizes and corrects for batch/capture effects | Pseudobulk counts | Batch-corrected expression matrix |

## Notes

Scripts 01 and 02 were run on a high-performance computing cluster. Update the file paths at the top of each script before running.
