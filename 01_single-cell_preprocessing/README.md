# 01. Single-cell preprocessing

Scripts to process raw scRNA-seq reads into a batch-corrected pseudobulk gene expression matrix. Run the scripts in numbered order.

| Script | Description | Input | Output |
|--------|-------------|-------|--------|
| `01_alignment.sh` | Aligns reads and generates count matrices using cellranger | Raw FASTQ files, config files | Gene-by-cell count matrices |
| `02_demultiplexing.sh` | Assigns cells to individuals using genotype data with demuxalot | Count matrices, genotype VCF | Cell-to-individual assignments |
| `03_cell_QC_and_annotation.R`, `03_cell_QC_and_annotation.sh` | Removes low-quality cells and doublets, assigns cell type labels from reference dataset | Count matrices, reference dataset | Filtered and annotated cells |
| `04_merge.R`, `04_merge.sh` | Creates a single Seurat file containing all captures | Individual Seurat objects by capture | Merged Seurat object |
| `05_celltype_annotation.R`, `05_celltype_annotation.sh` | Refines cell type annotations | Merged Seurat object | Finalized list of high-quality cells and their annotations |
| `06_integration_psuedobulk.R`, `06_integration_psuedobulk.sh` | Performs final integration and creates pseudobulk matrices | Finalized list of high-quality cells, merged Seurat object | UMAP coordinates, pseudobulk counts matrices |
| `07_capture_corrected_expression_matrices.R` | Corrects for technical capture effects in gene expression matrix | Pseudobulk counts matrices | Capture-corrected counts matrices |

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

Cells in each capture were assigned to individual donors with Demuxalot (v0.4.1), run through the Demuxafy Singularity container (v3.0.0), using donor genotypes generated as described in the Methods. Run this step after alignment has finished for a capture.

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

Before running, edit the paths in the "User settings" section at the top of `03_cell_QC_and_annotation.R`, as well as the `--partition` and `--account` lines in `03_cell_QC_and_annotation.sh`, to match your system.

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

Before running, edit the paths in the "User settings" section at the top of `04_merge.R`, as well as the `--partition` and `--account` lines in `04_merge.sh`, to match your system.

To merge into single Seurat object:

```
sbatch 04_merge.sh
```

## `05_celltype_annotation.R`

This script refines the cell type assignments for each of the eight broad cell types (B, CD4T, CD8T, DC, monocyte, NK, otherT, other). For each cell type, cells are normalized, integrated across conditions with RPCA, and clustered. Cells are kept if their reference cell type labels agree at all three resolutions (celltype.l1, l2, and l3) and they belong to a cluster retained after manual inspection of marker genes.

Before running, edit the paths in the "User settings" section at the top of `05_celltype_annotation.R`, as well as the `--partition` and `--account` lines in `05_celltype_annotation.sh`, to match your system.

To run:

```
mkdir -p logs
for ct in B CD4T CD8T DC monocyte NK otherT other; do
  sbatch --job-name=annot_${ct} 05_celltype_annotation.sh ${ct}
done
```

## `06_integration_pseudobulk.R`

Integrates all retained cells together (1,396,931 cells) and computes UMAP coordinates. It then sums counts per donor and condition (pseudobulk) for the five main cell types: CD4T, CD8T, NK, monocyte, and B. The influenza A gene "NA" is read into R as a missing value and dropped when Seurat objects are created, so its counts are read directly from the cellranger output and added to the pseudobulk count matrices as `IAV_NA`. A finalized list of retained cells and their UMAP coordinates are written to `final_pbmc_seurat_south_american_innate_metadata.rds` and the raw pseudobulk counts matrices are written to `raw_counts_list_unfiltered.rds`, both of which are available on Zenodo (https://doi.org/10.5281/zenodo.XXXXXXX).

Before running, edit the paths in the "User settings" section at the top of `06_integration_pseudobulk.R`, as well as the `--partition` and `--account` lines in `06_integration_pseudobulk.sh`, to match your system.

To run:

```
sbatch 06_pbmc_integration_pseudobulk.sh
```

## `07_capture_corrected_expression_matrices.R`



## Notes

Scripts 01-06 were run on a high-performance computing cluster. Update the file paths at the top of each script before running.
