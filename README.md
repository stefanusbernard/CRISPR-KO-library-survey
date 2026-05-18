# CRISPR-KO Library Survey

Analysis of sgRNA alignment quality and gene coverage across five human genome-wide CRISPR-KO libraries (Avana, Brunello, TKOv3, Yusa, and Jacquere), aligned against the T2T-CHM13v2.0 genome assembly.

---

## Repository structure

```
CRISPR-KO-library-survey/
│
├── scripts/
│   ├── 00_CRISPR-KO-library-survey-functions.R     # Shared R functions
│   ├── Fig02_AB_Fig03A_library_overview.Rmd         # → Fig2A, Fig2B, Fig3A
│   ├── Fig02_CJ_lfc_stratification.Rmd              # → Fig2C–D, SuppFig2A–J
│   ├── Fig03_BE_cds_pam_analysis.Rmd                # → Fig3B–E, SuppFig3A–D
│   ├── Fig04_AB_rescuing_guides.Rmd                 # → Fig4A–B
│   ├── LibrarySurvey_install_req_packages.R         # Package installer
│   └── data_preparation_scripts/
│       ├── sgrna_off_target_classification.Rmd      # Classifies sgRNA alignments
│       ├── normalize_lfc_utils.py                   # Shared normalization functions
│       └── normalize_{library}_lfc.ipynb            # One notebook per library
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
│   ├── read_count_data/             # Raw read counts per library/screen
│   ├── removed_genes_survey/        # Spacer/CDS/PAM data for removed-gene analysis
│   └── T2T_data/                    # T2T-CHM13v2.0 TxDb SQLite database
│
├── figure/                          # Output figures (.png)
└── results/                         # Supplementary tables (.xlsx)
```

---

## Running the analysis

### Prerequisites

Install R packages:

```r
source("scripts/LibrarySurvey_install_req_packages.R")
```

Install Python dependencies (for LFC normalization):

```bash
conda install pandas numpy matplotlib seaborn jupyter
```

### Step 1 — Classify sgRNA alignments

Run `scripts/data_preparation_scripts/sgrna_off_target_classification.Rmd` once per library.

**Requires:** alignment CSVs from GuideRefine (`data/sgrna_lfc_data/{library}_data/*_aln.csv`) and disposed sgRNA lists (`data/guiderefine_output/T2T-CHM13/*_disposed_sgRNAs.tsv`).

**Produces:** `data/sgrna_lfc_data/{library}_data/*_sgrna_alignment_classification_T2T.csv`

### Step 2 — Normalize LFC

Run each `scripts/data_preparation_scripts/normalize_{library}_lfc.ipynb` in Jupyter.

**Requires:** raw LFC files in `data/sgrna_lfc_data/{library}_data/`.

**Produces:** `data/sgrna_lfc_data/output_normalized/*_normalized_lfc.csv`

### Step 3 — Generate figures

Knit each `.Rmd` in `scripts/` from RStudio or via `rmarkdown::render()`.

| Script | Output |
|---|---|
| `Fig02_AB_Fig03A_library_overview.Rmd` | Fig2A, Fig2B, Fig3A |
| `Fig02_CJ_lfc_stratification.Rmd` | Fig2C–D, SuppFig2A–J |
| `Fig03_BE_cds_pam_analysis.Rmd` | Fig3B–E, SuppFig3A–D |
| `Fig04_AB_rescuing_guides.Rmd` | Fig4A–B |

---

## External dependency

Scripts source shared functions from the **GuideRefine** project:

```r
source("../../GuideRefine/GuideRefine_functions.R")
```

Update line 16 of `scripts/00_CRISPR-KO-library-survey-functions.R` if GuideRefine is at a different path.

---

## Data sources

See [`data/README.md`](data/README.md) for LFC file provenance, and [`data/annotation/README.md`](data/annotation/README.md) for genome annotation file sources.

> Large files (`.csv`, `.tsv`, `.xlsx`, `.gz`) are git-ignored. Re-download or copy from GuideRefine if cloning on a new machine.
