# A survey of multi-targeting and off-targeting sgRNAs across CRISPR-KO libraries

Analysis of sgRNA alignment quality and gene coverage across five human genome-wide CRISPR-KO libraries (Avana, Brunello, TKOv3, Yusa/Project Score, and Jacquere), aligned against the T2T-CHM13v2.0 genome assembly.

---

## Repository structure

```
CRISPR-KO-library-survey/
│
├── scripts/
│   ├── 00_CRISPR-KO-library-survey-functions.R   # Shared R functions (sourced by all figure scripts)
│   │
│   ├── Fig02_AB_Fig03A_library_overview.Rmd       # → Fig2A, Fig2B, Fig3A
│   ├── Fig02_CJ_lfc_stratification.Rmd            # → Fig2C–D, SuppFig2A–J
│   ├── Fig03_BE_cds_pam_analysis.Rmd              # → Fig3B–E, SuppFig3A–D
│   │
│   └── data_prep/                                 # Upstream scripts — produce data/, not figures
│       ├── classify_sgrna.Rmd                     # Classifies sgRNA alignments per library
│       ├── normalize_lfc_utils.py                 # Shared Python normalization functions
│       ├── normalize_avana_lfc.ipynb
│       ├── normalize_brunello_lfc.ipynb
│       ├── normalize_tkov3_lfc.ipynb
│       ├── normalize_yusa_lfc.ipynb
│       └── normalize_jacquere_lfc.ipynb
│
├── data/
│   ├── library_data/                              # Original sgRNA library files (.tsv)
│   ├── sgrna_lfc_data/                            # Raw LFC data, alignment CSVs, normalized output
│   ├── removed_genes_survey/                      # Spacer/CDS/PAM data for removed-gene analysis
│   ├── T2T_data/                                  # T2T-CHM13v2.0 TxDb SQLite database
│   └── biogrid_crispr_ko_data/                    # BioGRID ORCS CRISPR screen data
│
├── figure/
│   ├── Fig02_analysis_off-target-sgrna/           # All Fig2 and SuppFig2 outputs
│   └── Fig03_analysis_removed_genes/              # All Fig3 and SuppFig3 outputs
│
└── results/                                       # Supplementary tables (.xlsx)
```

---

## Analysis pipeline

The pipeline runs in three stages. Each stage depends on the outputs of the previous one.

### Stage 0 — sgRNA alignment (external project)

Alignment of sgRNA spacer sequences against T2T-CHM13v2.0 is performed by the **GuideRefine** project (located at `../../GuideRefine/` relative to this repository). GuideRefine produces per-library alignment files at:

```
../../GuideRefine/object_intermediate/T2T-CHM13_bsgenome/*_aln.csv
../../GuideRefine/output_cleaning/T2T-CHM13/*_disposed_sgRNAs.tsv
../../GuideRefine/output_cleaning/T2T-CHM13/*_full_report.xlsx
```

These files are required by `data_prep/classify_sgrna.Rmd` and `Fig02_AB_Fig03A_library_overview.Rmd`.

---

### Stage 1 — Data preparation (`scripts/data_prep/`)

Run these once to generate the processed data files consumed by figure scripts.

#### 1a. Classify sgRNA alignments

**Script:** `data_prep/classify_sgrna.Rmd`

Reads GuideRefine alignment output and classifies each sgRNA as: perfect on-target, single mismatch, PAM-distal mismatch, multi-target, or non-targeting.

**Outputs** → `data/sgrna_lfc_data/{library}_data/*_sgrna_alignment_classification_T2T.csv`

#### 1b. Normalize sgRNA Log2 fold-change (LFC)

**Scripts:** `data_prep/normalize_{library}_lfc.ipynb` (run in Jupyter; shared functions in `normalize_lfc_utils.py`)

Applies median-of-medians normalization followed by scaling to the absolute mean LFC of Hart (2014) common essential genes.

**Outputs** → `data/sgrna_lfc_data/output_normalized/*_normalized_lfc.csv`

---

### Stage 2 — Figure generation (`scripts/`)

Each script is self-contained and sources `00_CRISPR-KO-library-survey-functions.R`. Run by knitting the `.Rmd` file in RStudio, or via `rmarkdown::render()`.

| Script | Reads | Produces |
|---|---|---|
| `Fig02_AB_Fig03A_library_overview.Rmd` | alignment classifications, GuideRefine full reports | Fig2A, Fig2B, Fig3A; `removed_genes_all_library.csv` |
| `Fig02_CJ_lfc_stratification.Rmd` | alignment classifications, normalized LFC | Fig2C–D, SuppFig2A–J |
| `Fig03_BE_cds_pam_analysis.Rmd` | `removed_genes_all_library.csv`, spacer/CDS data, T2T TxDb | Fig3B–E, SuppFig3A–D; `supplementary_table_3_520_genes.xlsx` |

**Note:** `Fig03_BE_cds_pam_analysis.Rmd` contains a one-time setup chunk (labelled `create-txdb`, `eval=FALSE`) that builds the T2T-CHM13v2.0 TxDb SQLite database from the GFF annotation file. This only needs to be run once; the database is then saved to `data/T2T_data/TxDb_T2T-CHM13v2.sqlite` and loaded automatically on subsequent runs.

---

## Figure index

### Figure 2 — Off-target sgRNA analysis (`figure/Fig02_analysis_off-target-sgrna/`)

| File | Description |
|---|---|
| `Fig2A_percentage_off_target_T2T.png` | Percentage of multi-target, single-mismatch, and PAM-distal mismatch sgRNAs per library |
| `Fig2B_percentage_removed_sgrna_T2T.png` | Percentage of sgRNAs removed after library refinement |
| `Fig2C_combined_library_multi_target.png` | LFC distribution: perfect vs multi-target sgRNAs (all libraries) |
| `Fig2D_combined_library_single_mismatch.png` | LFC distribution: perfect vs single-mismatch sgRNAs (all libraries) |
| `SuppFig2A–E_*_multi_target.png` | Per-library LFC stratified by number of perfect on-target alignments |
| `SuppFig2F–J_*_single_mismatch.png` | Per-library LFC stratified by number of single-mismatch alignments |

### Figure 3 — Removed genes analysis (`figure/Fig03_analysis_removed_genes/`)

| File | Description |
|---|---|
| `Fig3A_library_genes_bar_chart_T2T.png` | Percentage of protein-coding genes with < 3 sgRNAs per library |
| `Fig3B_removed_genes_all_libraries.png` | UpSet plot of gene overlap across libraries |
| `Fig3C_log10_cds_length_removed_genes.png` | CDS length: removed vs control genes (log10 boxplot) |
| `Fig3D_log10_pam_removed_genes.png` | PAM-site count: removed vs control genes (log10 boxplot) |
| `Fig3E_pam_density_distribution.png` | PAM-site density (per kb CDS): removed vs control genes |
| `SuppFig3A–B_cds_length_*.png` | CDS length distribution barplots |
| `SuppFig3C–D_spacer_pam_*.png` | Protospacer count distribution barplots |

---

## Prerequisites

### R packages

```r
install.packages(c("tidyverse", "ggbreak", "ggsignif", "ggpubr", "rstatix",
                   "clinfun", "effsize", "UpSetR", "HGNChelper", "openxlsx",
                   "RColorBrewer", "paletteer", "VennDiagram"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "txdbmaker", "crisprVerse",
                       "BSgenome.Hsapiens.NCBI.T2TCHM13v2.0",
                       "BSgenome.Hsapiens.UCSC.hg38",
                       "BSgenome", "GenomeInfoDbData"))
```

### Python (for LFC normalization notebooks)

Python ≥ 3.9 with `pandas`, `numpy`, `matplotlib`, `seaborn`. Recommended via conda:

```bash
conda install pandas numpy matplotlib seaborn jupyter
```

---

## Data sources

See [`data/README.md`](data/README.md) for full provenance of all raw data files (LFC matrices, library TSVs, BioGRID ORCS data).
