# Data

Data files are not stored in this repository. There are two ways to get them:

1. **Start from the processed data (recommended).** Download `data.zip` from Zenodo (https://doi.org/10.5281/zenodo.23425836) into the main repository folder and unzip it. You can then run the downstream analyses directly without the preprocessing steps.
2. **Reproduce the preprocessing.** Download the raw data from GEO (GSE351087), and run the scripts in `01_single-cell_preprocessing/`. They write their outputs to the same folders.

| Folder | Contents | Source |
|--------|----------|--------|
| `metadata/` | Sample metadata | Zenodo |
| `integration/` | List of retained cells after QC, raw counts matrices | Zenodo or Scripts 05–06 |
| `counts_and_weights/` | Capture-corrected expression matrices and voom weights | Zenodo or Script 07 |
