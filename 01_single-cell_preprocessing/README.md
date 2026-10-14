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

## `01_alignment.sh`

Reads from each of the 63 captures were aligned separately with cellranger multi (v9.0.1) to a combined reference of the human genome (hg38), influenza A virus, and *Trypanosoma cruzi*. The settings for each capture, including the sequencing libraries and lanes used, are in `configs/<capture>_config.csv`. The full list of captures is in `captures.txt`.

Before running, edit the paths in the "User settings" section at the top of `01_alignment.sh`, as well as the `--partition` and `--account` lines, to match your system. The config files use the placeholders `REFERENCE_PATH` and `FASTQ_DIR`, which the script fills in automatically from those settings.

To align a single capture:

```
mkdir -p logs
sbatch --job-name=align_IAV-b1 01_alignment.sh IAV-b1
```

To align all captures:

```
mkdir -p logs
for capture in $(cat captures.txt); do
  sbatch --job-name=align_${capture} 01_alignment.sh ${capture}
done
```

Log files for each capture are written to `logs/`. Note that `--jobmode=slurm` requires cellranger's SLURM job template to be configured on your cluster.

## Notes

Scripts 01 and 02 were run on a high-performance computing cluster. Update the file paths at the top of each script before running.
