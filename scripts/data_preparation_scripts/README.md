# Data Preparation Scripts

These scripts generate the intermediate data files consumed by the main figure scripts. Run them in order before knitting any `Fig*.Rmd` or `SuppTable*.Rmd`.

---

## Execution order

```
(only when adding a new library)
Step -1 → reformat_sgRNA_library.rmd or reformat_Jacquere_library.Rmd

Step 0 → build_restricted_library.Rmd          (once)
          → then run GuideRefine (external) on each restricted library
Step 1 → sgrna_off_target_classification.Rmd   (once per library)
Step 2 → normalize_{library}_lfc.ipynb         (once per library, in
          ../sgrna_log_fold_change_scripts/ — see that folder's README)
```

---

## Scripts

### `reformat_sgRNA_library.rmd` / `reformat_Jacquere_library.Rmd` (R)

Convert a raw, publicly-downloaded sgRNA library into GuideRefine's 3-column format (`sgRNA`, `spacer`, `gene`, no header). `reformat_sgRNA_library.rmd` handles standard single-gene-per-row libraries (Avana, Brunello, TKOv3, Yusa); `reformat_Jacquere_library.Rmd` additionally resolves Jacquere's multi-gene `GENE1|GENE2` fields. Only needed when adding a library — the 5 libraries used in this study are already reformatted.

**Requires:** `public_crispr_library/raw/*`

**Produces:** `public_crispr_library/processed/*.tsv`

See [`public_crispr_library/README.md`](../../public_crispr_library/README.md) for details.

---

### `build_restricted_library.Rmd` (R)

Splits out control sgRNAs and restricts targeting guides to genes that are both (a) HGNC protein-coding and (b) have a RefSeq Select-tagged canonical transcript in T2T-CHM13.

**Requires:** `public_crispr_library/processed/*.tsv`, an HGNC complete-set snapshot (`data/hgnc_complete_set_*.txt`), and `data/annotation/T2T-CHM13v2.0_gene_annot_granges.rds`.

**Produces:** `public_crispr_library/restricted_library/*.tsv` — this is the actual GuideRefine input used in this study.

---

### `sgrna_off_target_classification.Rmd` (R)

Classifies every sgRNA alignment in a library into one of five categories:
`multi-target guides`, `single mismatch`, `pam-distal single mismatch`, `pam-distal double mismatch`, or `perfect`.

**Requires:**
- Alignment CSV from GuideRefine: `data/sgrna_lfc_data/{library}_data/*_aln.csv`
- Disposed sgRNA list from GuideRefine: `data/guiderefine_output/Jul2026_T2T-CHM13/*_disposed_sgRNAs.tsv`
- T2T-CHM13v2.0 annotation object: `data/annotation/T2T-CHM13v2.0_gene_annot_granges.rds`
- hg38 CCDS file (hg38 mode only): `data/annotation/CCDS.20221027.txt`

**Produces:** `data/sgrna_lfc_data/{library}_data/*_sgrna_alignment_classification_T2T.csv`

Run once per library (avana, brunello, tkov3, yusa, jacquere). The script supports both T2T-CHM13 and hg38 — two separate code chunks are provided; run the relevant one and skip the other.

---

LFC normalization (`normalize_{library}_lfc.ipynb`, `normalize_lfc_utils.py`, the `depmap_ppi` conda environment) now lives in [`../sgrna_log_fold_change_scripts/`](../sgrna_log_fold_change_scripts/README.md).
