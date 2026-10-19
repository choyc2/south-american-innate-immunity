# Data

Data files are not stored in this repository. Download `data.zip` from Zenodo (https://doi.org/10.5281/zenodo.23425836) into the main repository folder and unzip it. This creates three folders:

| Folder | Contents | Source |
|--------|----------|--------|
| `metadata/` | Sample metadata | Zenodo |
| `integration/` | List of retained cells after QC and their UMAP coordinates, raw counts matrices | Zenodo or Scripts 05–06 |
| `counts_and_weights/` | Capture-corrected expression matrices and voom weights | Zenodo or Script 07 |

With these files you can:

1. **Start from the processed data (recommended).** All files required to run the downstream main analyses are provided in `data.zip` so that you do not have to run the preprocessing steps yourself.
2. **Reproduce the preprocessing.** Download the raw scRNA-seq data from GEO (GSE351087), and run the scripts in `01_single-cell_preprocessing/`. Note that you still need to download `data.zip` since preprocessing requires the files in `metadata/`. The scripts will regenerate the files in `integration/` and `counts_and_weights/`, overwriting the downloaded copies.
