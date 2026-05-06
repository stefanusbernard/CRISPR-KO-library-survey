# A survey of multi-targeting and off-targeting sgRNAs across five genome-wide CRISPR-Cas9 knockout screen’s sgRNA libraries 

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
│   ├── annotation/                                # Genome annotation files (T2T GRanges RDS, CCDS, GFF)
│   ├── library_data/                              # Original sgRNA library files (.tsv)
│   ├── sgrna_lfc_data/                            # Raw LFC, alignment CSVs, normalized output
│   │   ├── avana_data/
│   │   ├── brunello_data/
│   │   ├── tkov3_data/
│   │   ├── yusa_data/
│   │   ├── jacquere_data/
│   │   ├── hg38_alignment/                        # hg38 *_aln.csv files (git-ignored, local only)
│   │   └── output_normalized/
│   ├── guiderefine_output/                        # Full reports and disposed sgRNA lists from GuideRefine
│   │   ├── T2T-CHM13/
│   │   │   └── old_result_pam_distal_mismatch_only/
│   │   └── hg38/
│   │       └── old_public_library_result/
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

Alignment of sgRNA spacer sequences against T2T-CHM13v2.0 is performed by the **GuideRefine** project (located at `../../GuideRefine/` relative to this repository).

The output files produced by GuideRefine have been copied into this repository under `data/` and all scripts now read from there directly. The only remaining dependency on the GuideRefine project is its shared R function library:

```r
source("../../GuideRefine/GuideRefine_functions.R")  # sourced by 00_CRISPR-KO-library-survey-functions.R
```

This library is intentionally kept in GuideRefine to avoid maintaining two copies. If GuideRefine is not available at `../../GuideRefine/`, update the `source()` path on line 16 of `scripts/00_CRISPR-KO-library-survey-functions.R`.

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
| `Fig02_AB_Fig03A_library_overview.Rmd` | alignment classifications, `data/guiderefine_output/` full reports | Fig2A, Fig2B, Fig3A; `removed_genes_all_library.csv` |
| `Fig02_CJ_lfc_stratification.Rmd` | alignment classifications, normalized LFC | Fig2C–D, SuppFig2A–J |
| `Fig03_BE_cds_pam_analysis.Rmd` | `removed_genes_all_library.csv`, spacer/CDS data, T2T TxDb | Fig3B–E, SuppFig3A–D; `supplementary_table_3_520_genes.xlsx` |

**Note:** `Fig03_BE_cds_pam_analysis.Rmd` contains a one-time setup chunk (labelled `create-txdb`, `eval=FALSE`) that builds the T2T-CHM13v2.0 TxDb SQLite database from the GFF annotation file. This only needs to be run once; the database is then saved to `data/T2T_data/TxDb_T2T-CHM13v2.sqlite` and loaded automatically on subsequent runs.

---

## GuideRefine data (bundled locally)

All files originally produced by the GuideRefine project have been copied into `data/` so scripts run without requiring the external repository. The table below shows where each file now lives and its git-tracking status.

> **Note:** Large CSV and TSV files are excluded from git by `.gitignore` (`*.csv`, `*.tsv`, `*.xlsx`, `*.gz`). They are local copies only — copy them from GuideRefine again if you clone this repository on a new machine.

### Genome annotation — `data/annotation/`

| File | Size | Git-tracked | Used by |
|---|---|---|---|
| `T2T-CHM13v2.0_gene_annot_granges.rds` | 1.4 MB | ✅ yes | `data_prep/classify_sgrna.Rmd` |
| `CCDS.20221027.txt` | 10 MB | ❌ no | `data_prep/classify_sgrna.Rmd` (hg38 chunk) |
| `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz` | 76 MB | ❌ no | `Fig03_BE_cds_pam_analysis.Rmd` (one-time `create-txdb` chunk, `eval=FALSE`) |

### sgRNA library sequences — `data/library_data/original_library/`

Five library TSV files (Avana, Brunello, TKOv3, Yusa, Jacquere). Git-ignored. Used by `data_prep/classify_sgrna.Rmd` and `data_prep/normalize_*_lfc.ipynb`.

### T2T-CHM13v2.0 alignment CSVs — `data/sgrna_lfc_data/{library}_data/`

| File | Destination folder | Size |
|---|---|---|
| `avana_library_aln.csv` | `avana_data/` | 33 MB |
| `broadgpp-brunello-library-contents_aln.csv` | `brunello_data/` | 110 MB |
| `tkov3_guide_sequence_aln.csv` | `tkov3_data/` | 11 MB |
| `yusa_hcrispr_ko_grnas_aln.csv` | `yusa_data/` | 22 MB |
| `Jacquere_PerGuideAnnotations_Quota4_aln.csv` | `jacquere_data/` | 18 MB |

All git-ignored. Used by `data_prep/classify_sgrna.Rmd`.

### hg38 alignment CSVs — `data/sgrna_lfc_data/hg38_alignment/`

Same five filenames as above (35–118 MB each). Git-ignored. Used only by the hg38 chunk of `data_prep/classify_sgrna.Rmd`.

### Disposed sgRNA lists (T2T-CHM13v2.0) — `data/guiderefine_output/T2T-CHM13/old_result_pam_distal_mismatch_only/`

Five `*_disposed_sgRNAs.tsv` files (one per library). Git-ignored. Used by `data_prep/classify_sgrna.Rmd`.

### Disposed sgRNA lists (hg38) — `data/guiderefine_output/hg38/old_public_library_result/`

Same five filenames. Git-ignored. Used only by the hg38 chunk of `data_prep/classify_sgrna.Rmd`.

### Library refinement reports — `data/guiderefine_output/T2T-CHM13/`

Five `*_full_report.xlsx` files (one per library). Git-ignored. Used by `Fig02_AB_Fig03A_library_overview.Rmd`.

### Remaining GuideRefine dependency

The only file **not** copied is `GuideRefine_functions.R`. It is sourced directly from the GuideRefine repository to avoid maintaining two copies. If GuideRefine is unavailable, update line 16 of `scripts/00_CRISPR-KO-library-survey-functions.R`:

```r
source("D:/GitHub/GuideRefine/GuideRefine_functions.R")  # update this path
```

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
