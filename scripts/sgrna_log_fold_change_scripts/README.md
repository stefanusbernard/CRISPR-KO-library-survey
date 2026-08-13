# sgRNA Log-Fold-Change Scripts

Normalizes raw per-guide log-fold-change (LFC) data for each library. Run after `sgrna_off_target_classification.Rmd` (see [`../data_preparation_scripts/README.md`](../data_preparation_scripts/README.md) for the earlier pipeline steps).

---

## Scripts

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

**Requires:** `data/sgrna_lfc_data/avana_data/` (raw LFC files), `public_crispr_library/restricted_library/avana_library.tsv`

**Produces:** `data/sgrna_lfc_data/output_normalized/avana_normalized_lfc.csv`

---

### `normalize_brunello_lfc.ipynb` (Python)

Normalizes raw LFC data for the **Brunello** library (A375 cell line, Sanson et al. 2018).

**Requires:** `data/sgrna_lfc_data/brunello_data/` (raw LFC files), `public_crispr_library/restricted_library/broadgpp-brunello-library-contents.tsv`

**Produces:** `data/sgrna_lfc_data/output_normalized/brunello_a375_normalized_lfc.csv`

---

### `normalize_jacquere_lfc.ipynb` (Python)

Normalizes raw LFC data for the **Jacquere** library.

**Requires:** `data/sgrna_lfc_data/jacquere_data/` (raw LFC files), `public_crispr_library/restricted_library/Jacquere_PerGuideAnnotations_Quota4.tsv`

**Produces:** `data/sgrna_lfc_data/output_normalized/jacquere_normalized_lfc.csv`

---

### `normalize_tkov3_lfc.ipynb` (Python)

Normalizes raw LFC data for the **TKOv3** library (RPE-1 cell line).

**Requires:** `data/sgrna_lfc_data/tkov3_data/` (raw LFC files), `public_crispr_library/restricted_library/tkov3_guide_sequence.tsv`

**Produces:** `data/sgrna_lfc_data/output_normalized/tkov3_rpe1_normalized_lfc.csv`

---

### `normalize_yusa_lfc.ipynb` (Python)

Normalizes raw LFC data for the **Yusa (Project Score)** library.

**Requires:** `data/sgrna_lfc_data/yusa_data/` (raw LFC files), `public_crispr_library/restricted_library/yusa_hcrispr_ko_grnas.tsv`

**Produces:** `data/sgrna_lfc_data/output_normalized/yusa_normalized_lfc.csv`

---

## Normalization method

All five notebooks apply **essential-gene scaling** via `normalize_lfc_utils.py`'s `scale_essential()` — LFC values are divided by the mean absolute LFC of Hart 2014 core-essential genes (`data/sgrna_lfc_data/constitutive_core_essential_hart_2014.csv`), anchoring the scale so that essential-gene depletion corresponds to approximately −1.

Only **Brunello, TKOv3, and Jacquere** chain a median-MAD normalization step (`normalize_column()` — each LFC column centred by its median and scaled by its MAD) before essential-gene scaling. **Avana and Yusa** compute `normalize_column()` too, but only for the before/after sanity-check plot (`sanity_check_scale_essential()`) — their actual output comes from essential-gene scaling applied directly to the raw per-column LFC, skipping the median-MAD step. This is intentional (Avana and Yusa are both large multi-cell-line DepMap-format matrices, treated identically), not a bug — but it means the two library groups aren't normalized the same way, worth knowing if comparing raw scale/spread across libraries.

The notebooks were run using the `depmap_ppi` conda environment (Python 3.12.7, Anaconda distribution). The full environment is exported to [`depmap_ppi_environment.yml`](depmap_ppi_environment.yml) in this directory — recreate it with:

```bash
conda env create -f scripts/sgrna_log_fold_change_scripts/depmap_ppi_environment.yml
conda activate depmap_ppi
```

Key package versions: pandas 2.2.2, numpy 1.26.4, matplotlib 3.9.2, seaborn 0.13.2, jupyterlab 4.2.5.
