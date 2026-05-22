# Data Preparation Scripts

These scripts generate the intermediate data files consumed by the main figure scripts. Run them in order before knitting any `Fig*.Rmd` or `SuppTable*.Rmd`.

---

## Execution order

```
Step 1 → sgrna_off_target_classification.Rmd   (once per library)
Step 2 → normalize_{library}_lfc.ipynb         (once per library)
```

---

## Scripts

### `sgrna_off_target_classification.Rmd` (R)

Classifies every sgRNA alignment in a library into one of five categories:
`multi-target guides`, `single mismatch`, `pam-distal single mismatch`, `pam-distal double mismatch`, or `perfect`.

**Requires:**
- Alignment CSV from GuideRefine: `data/sgrna_lfc_data/{library}_data/*_aln.csv`
- Disposed sgRNA list from GuideRefine: `data/guiderefine_output/May2026_T2T-CHM13/*_disposed_sgRNAs.tsv`
- T2T-CHM13v2.0 annotation object: `data/annotation/T2T-CHM13v2.0_gene_annot_granges.rds`
- hg38 CCDS file (hg38 mode only): `data/annotation/CCDS.20221027.txt`

**Produces:** `data/sgrna_lfc_data/{library}_data/*_sgrna_alignment_classification_T2T.csv`

Run once per library (avana, brunello, tkov3, yusa, jacquere). The script supports both T2T-CHM13 and hg38 — two separate code chunks are provided; run the relevant one and skip the other.

---

### `normalize_lfc_utils.py` (Python — shared utility)

Shared helper functions imported by all five normalization notebooks. Not run directly.

Provides:
- `normalize_column()` — median-centres and MAD-scales a single LFC column.
- `scale_essential()` — scales LFC values so that the mean absolute LFC of Hart 2014 core-essential genes equals 1.
- `sanity_check_scale_essential()` — boxplot QC showing LFC distributions before and after scaling.
- `mean_row()` — averages LFC across cell lines per sgRNA.

---

### `normalize_avana_lfc.ipynb` (Python)

Normalizes raw LFC data for the **Avana** library.

**Requires:** `data/sgrna_lfc_data/avana_data/` (raw LFC files)

**Produces:** `data/sgrna_lfc_data/output_normalized/avana_normalized_lfc.csv`

---

### `normalize_brunello_lfc.ipynb` (Python)

Normalizes raw LFC data for the **Brunello** library (A375 cell line, Sanson et al. 2018).

**Requires:** `data/sgrna_lfc_data/brunello_data/` (raw LFC files)

**Produces:** `data/sgrna_lfc_data/output_normalized/brunello_a375_normalized_lfc.csv`

---

### `normalize_jacquere_lfc.ipynb` (Python)

Normalizes raw LFC data for the **Jacquere** library.

**Requires:** `data/sgrna_lfc_data/jacquere_data/` (raw LFC files)

**Produces:** `data/sgrna_lfc_data/output_normalized/jacquere_normalized_lfc.csv`

---

### `normalize_tkov3_lfc.ipynb` (Python)

Normalizes raw LFC data for the **TKOv3** library (RPE-1 cell line).

**Requires:** `data/sgrna_lfc_data/tkov3_data/` (raw LFC files)

**Produces:** `data/sgrna_lfc_data/output_normalized/tkov3_rpe1_normalized_lfc.csv`

---

### `normalize_yusa_lfc.ipynb` (Python)

Normalizes raw LFC data for the **Yusa (Project Score)** library.

**Requires:** `data/sgrna_lfc_data/yusa_data/` (raw LFC files)

**Produces:** `data/sgrna_lfc_data/output_normalized/yusa_normalized_lfc.csv`

---

## Normalization method (all notebooks)

Each notebook applies the same two-step normalization via `normalize_lfc_utils.py`:

1. **Median-MAD normalization** — each LFC column is centred by its median and scaled by its MAD, making distributions comparable across screens.
2. **Essential-gene scaling** — LFC values are divided by the mean absolute LFC of Hart 2014 core-essential genes (`data/sgrna_lfc_data/constitutive_core_essential_hart_2014.csv`), anchoring the scale so that essential-gene depletion corresponds to approximately −1.

The notebooks were run using the `depmap_ppi` conda environment (Python 3.12.7).
