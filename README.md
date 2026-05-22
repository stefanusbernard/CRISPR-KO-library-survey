# A survey of multi-targeting and off-targeting sgRNAs across five genome-wide CRISPR-Cas9 knockout screen's sgRNA libraries

Analysis of sgRNA alignment quality and gene coverage across five human genome-wide CRISPR-KO libraries (Avana, Brunello, TKOv3, Yusa, and Jacquere), aligned against the T2T-CHM13v2.0 genome assembly.

---

## Repository structure

```
CRISPR-KO-library-survey/
│
├── scripts/
│   ├── GuideRefine-functions.R                      # Local copy of GuideRefine shared functions
│   ├── CRISPR-KO-library-survey-functions.R         # Shared R functions for this project
│   ├── Fig02_AB_Fig03A_library_overview.Rmd         # → Fig2A, Fig2B, Fig3A
│   ├── Fig02_CJ_lfc_stratification.Rmd              # → Fig2C–D, SuppFig2A–J
│   ├── Fig03_BE_cds_pam_analysis.Rmd                # → Fig3B–E, SuppFig3A–D
│   ├── SuppTable03_ABC_aggregate_sgrna.Rmd          # → SuppTable3A–C, Fig4A–B
│   ├── LibrarySurvey_install_req_packages.R         # Package installer (run first)
│   ├── T2T_annotation_scripts/                      # One-time genome resource setup
│   │   ├── 01_T2T_gff_to_ccds_conversion.Rmd       # GFF → GRanges annotation object
│   │   └── 02_forge_T2T_bsgenome.Rmd               # Build BSgenome.Hsapiens.NCBI.T2TCHM13v2.0
│   └── data_preparation_scripts/
│       ├── sgrna_off_target_classification.Rmd      # Classifies sgRNA alignments
│       ├── normalize_lfc_utils.py                   # Shared normalization functions
│       ├── normalize_avana_lfc.ipynb
│       ├── normalize_brunello_lfc.ipynb
│       ├── normalize_jacquere_lfc.ipynb
│       ├── normalize_tkov3_lfc.ipynb
│       └── normalize_yusa_lfc.ipynb
│
├── data/
│   ├── annotation/                  # Genome annotation files (see annotation/README.md)
│   ├── library_data/
│   │   ├── original_library/        # Original sgRNA library TSV files
│   │   └── refined_library/         # Refined library TSV files (GuideRefine output)
│   ├── sgrna_lfc_data/              # Alignment CSVs and normalized LFC per library
│   │   ├── avana_data/
│   │   ├── brunello_data/
│   │   ├── tkov3_data/
│   │   ├── yusa_data/
│   │   ├── jacquere_data/
│   │   ├── hg38_alignment/
│   │   └── output_normalized/
│   ├── guiderefine_output/          # Full reports and disposed sgRNA lists
│   │   ├── T2T-CHM13/
│   │   ├── May2026_T2T-CHM13/
│   │   └── hg38/
│   ├── hits_change_data/            # DepMap gene effect and LFC data
│   └── read_count_data/             # Raw read counts per library/screen
│
└── figure/                          # Output figures (.png)
```

---

## Running the analysis

### Prerequisites

**Install R packages first** — run this before any other script:

```r
source("scripts/LibrarySurvey_install_req_packages.R")
```

This installs all required CRAN and Bioconductor packages. Run it once before proceeding to any of the steps below.

Install Python dependencies (for LFC normalization):

```bash
conda install pandas numpy matplotlib seaborn jupyter
```

### Step 0 — One-time genome setup (new machine only)

Run the scripts in `scripts/T2T_annotation_scripts/` once to build the genome resources. See [`scripts/T2T_annotation_scripts/README.md`](scripts/T2T_annotation_scripts/README.md) for download instructions and full details.

| Script | Produces |
|---|---|
| `01_T2T_gff_to_ccds_conversion.Rmd` | `data/annotation/T2T-CHM13v2.0_gene_annot_granges.rds` |
| `02_forge_T2T_bsgenome.Rmd` | `BSgenome.Hsapiens.NCBI.T2TCHM13v2.0` R package |

### Step 1 — Classify sgRNA alignments

Run `scripts/data_preparation_scripts/sgrna_off_target_classification.Rmd` once per library.

**Requires:** alignment CSVs from GuideRefine (`data/sgrna_lfc_data/{library}_data/*_aln.csv`) and disposed sgRNA lists (`data/guiderefine_output/T2T-CHM13/*_disposed_sgRNAs.tsv`).

**Produces:** `data/sgrna_lfc_data/{library}_data/*_sgrna_alignment_classification_T2T.csv`

### Step 2 — Normalize LFC

Run each `scripts/data_preparation_scripts/normalize_{library}_lfc.ipynb` in Jupyter.

**Requires:** raw LFC files in `data/sgrna_lfc_data/{library}_data/`.

**Produces:** `data/sgrna_lfc_data/output_normalized/*_normalized_lfc.csv`

### Step 3 — Generate figures and supplementary tables

Knit each `.Rmd` in `scripts/` from RStudio or via `rmarkdown::render()`.

| Script | Output |
|---|---|
| `Fig02_AB_Fig03A_library_overview.Rmd` | Fig2A, Fig2B, Fig3A |
| `Fig02_CJ_lfc_stratification.Rmd` | Fig2C–D, SuppFig2A–J |
| `Fig03_BE_cds_pam_analysis.Rmd` | Fig3B–E, SuppFig3A–D |
| `SuppTable03_ABC_aggregate_sgrna.Rmd` | SuppTable3A–C, Fig4A–B |

---

## Shared functions

All scripts source two function files at the project root:

```r
source(here::here("scripts/GuideRefine-functions.R"))
source(here::here("scripts", "CRISPR-KO-library-survey-functions.R"))
```

`scripts/GuideRefine-functions.R` is a local copy of the shared alignment and annotation functions from the GuideRefine pipeline. If you update GuideRefine, copy the updated file here to keep them in sync.

---

## Data sources

See [`data/README.md`](data/README.md) for LFC file provenance, and [`data/annotation/README.md`](data/annotation/README.md) for genome annotation file sources.

> Large files (`.csv`, `.tsv`, `.xlsx`, `.gz`) are git-ignored. Re-download or copy from GuideRefine if cloning on a new machine.

---

## Tools and versions

### R (4.5.2)

**CRAN packages**

| Package | Version | Purpose |
|---|---|---|
| tidyverse | 2.0.0 | Data wrangling and plotting (ggplot2, dplyr, readr, tidyr, purrr) |
| readxl | 1.5.0 | Read Excel files |
| openxlsx | 4.2.8.1 | Write Excel files |
| ggbreak | 0.1.7 | Axis breaks for ggplot2 |
| UpSetR | 1.4.0 | UpSet intersection plots |
| ggpubr | 0.6.3 | Publication-ready ggplot2 figures |
| ggsignif | 0.6.4 | Significance brackets for ggplot2 |
| HGNChelper | 0.8.15 | HGNC gene symbol validation and correction |
| here | 1.0.2 | Project-root-aware file paths |
| rmarkdown | 2.31 | R Markdown rendering |
| knitr | 1.51 | Dynamic report generation |
| pandoc | 0.2.0 | Document conversion (system tool via rmarkdown) |

**Bioconductor packages**

| Package | Version | Purpose |
|---|---|---|
| BiocManager | 1.30.27 | Bioconductor package installer |
| BSgenome | 1.78.0 | Infrastructure for full genome sequences |
| BSgenome.Hsapiens.UCSC.hg38 | 1.4.5 | hg38 genome sequence |
| BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 | 1.5.0 | T2T-CHM13v2.0 genome sequence |
| GenomicFeatures | 1.62.0 | Genomic annotation infrastructure |
| GenomeInfoDbData | 1.2.15 | Chromosome metadata |
| txdbmaker | 1.6.2 | Build TxDb objects from annotation sources |
| crisprVerse | 1.12.0 | CRISPR guide design framework |
| crisprBase | 1.14.0 | Core CRISPR data structures |

---

### Python (3.12.7)

Notebooks were run using the `depmap_ppi` conda environment (Python 3.12.7).

| Package | Purpose |
|---|---|
| pandas | Tabular data manipulation |
| numpy | Numerical operations (median, MAD normalisation) |
| matplotlib | Base plotting library |
| seaborn | Statistical visualisation |
| jupyter | Interactive notebook environment |
