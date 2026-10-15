# 01. Single-cell preprocessing

Scripts to process raw single-cell RNA-seq reads into a batch-corrected pseudobulk gene expression matrix. Run the scripts in numbered order.

| Script | Description | Input | Output |
|--------|-------------|-------|--------|
| `01_alignment.sh` | Aligns reads and generates count matrices using cellranger multi (10x Genomics, v9.0.1) | Raw FASTQ files, config files | Gene-by-cell count matrices |
| `02_demultiplexing.sh` | Assigns cells to individuals using genotype data with demuxalot (v0.4.1) | Count matrices, genotype VCF | Cell-to-individual assignments |
| `03_cell_QC_and_annotation.R`, `03_cell_QC_and_annotation.sh` | Removes low-quality cells and doublets, assigns cell type labels from reference dataset | Count matrices, reference dataset | Filtered cells |
| `04_merge.R`, `04_merge.sh` | Creates a single Seurat file containing all captures | Individuals capture Seurat objects | Merged Seurat object |
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

## `02_demultiplexing.sh`

Cells in each capture were assigned to individual donors with Demuxalot, run through the Demuxafy Singularity container (v3.0.0), using donor genotypes generated as described in the Methods. Run this step after alignment has finished for a capture.

Before running, edit the paths in the "User settings" section of `02_demultiplexing.sh`. `GENOTYPE_DIR` must contain, for each sequencing batch, a list of the donors in that batch (`batch<N>_inds.txt`) and their genotypes (`batch<N>.MAF0.05.GRCh38_v2.vcf.gz`). These files contain individual-level genetic data and are not included in this repository; see the [main README](https://github.com/choyc2/south-american-innate-immunity/blob/main/README.md) for data access. The batch is determined from the capture name (for example, `IAV-b1` uses batch 1).

To demultiplex a single capture:

```
mkdir -p logs
sbatch --job-name=demux_IAV-b1 02_demultiplexing.sh IAV-b1
```

To demultiplex all captures:

```
mkdir -p logs
for capture in $(cat captures.txt); do
  sbatch --job-name=demux_${capture} 02_demultiplexing.sh ${capture}
done
```

## `03_cell_QC_and_annotation.R`

This script creates a Seurat object for each capture, adds the Demuxalot donor assignments, and removes cells with ≤200 or ≥6,000 detected genes, ≥10% mitochondrial reads, or a doublet assignment. Cell types are then assigned at three resolutions by reference mapping to the PBMC atlas from [Hao et al. (2021)](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421005833%3Fshowall%3Dtrue), downloaded from [here](https://atlas.fredhutch.org/nygc/multimodal-pbmc/).

To run QC on a single capture:

```
mkdir -p logs
sbatch --job-name=qc_IAV-b1 03_cell_QC_and_annotation.sh IAV-b1

```

To run QC on all captures:

```
mkdir -p logs
for capture in $(cat captures.txt); do
  sbatch --job-name=qc_${capture} 03_cell_QC_and_annotation.sh ${capture}
done
```

## `04_merge.R`

Once each capture is converted to its own Seurat file, this script merges all captures into a single unnormalized Seurat object, adds condition, batch, and 10x chemistry to the metadata, and removes donors not included in this study.

To merge into single Seurat object:

```
sbatch 04_merge.sh
```

## Notes

Scripts 01-05 were run on a high-performance computing cluster. Update the file paths at the top of each script before running.
