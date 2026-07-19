# Supplementary Materials

This folder contains the supplementary data files companion for "A survey of multi-targeting and off-targeting sgRNAs across five genome-wide CRISPR-Cas9 knockout screen's sgRNA libraries"

---

## What is in this folder?

The files are organized into three groups:

1. **Supplementary materials 1–5** — tables describing how each guide RNA in each library aligns to the human genome
2. **Supplementary tables 1C and 1D** — summary tables about gene coverage across libraries
3. **GuideRefine outputs** — detailed quality filtering results for each library

---

## Supplementary Materials 1–5: sgRNA Alignment Classification

These five CSV (comma-separated) files each correspond to one of the five libraries evaluated in our study. For every guide RNA (sgRNA) in the library, the table records where it maps in the human genome and how that alignment was classified.

| File | Library |
|---|---|
| `supp_material_1_brunello_sgrna_alignment_classification_T2T.csv` | Brunello |
| `supp_material_2_tkov3_sgrna_alignment_classification_T2T.csv` | TKOv3 |
| `supp_material_3_yusa_sgrna_alignment_classification_T2T.csv` | Yusa |
| `supp_material_4_avana_sgrna_alignment_classification_T2T.csv` | Avana |
| `supp_material_5_jacquere_sgrna_alignment_classification_T2T.csv` | Jacquere |

Each file has the following columns:

| Column | What it means |
|---|---|
| `sgRNA` | A unique identifier for the guide RNA (e.g. `sgTP53_1`) |
| `spacer` | The 20-nucleotide sequence of the guide RNA |
| `gene` | The gene the guide was designed to target |
| `n_mismatches` | How many mismatches exist in the best genomic alignment (0 = perfect match) |
| `num_alignments` | Total number of places in the genome this guide aligns to |
| `alignment` | A summary classification of the alignment. Possible values: `perfect`, `single mismatch`, `pam-distal double mismatch`, `multi-target guides`, `non-targeting` |

---

## Supplementary Tables 1C and 1D

These Excel files summarize gene coverage across the five libraries.

**`supp_table_1c_aggregated_on_target_guides_represented_genes.xlsx`**
The mini-composite library derived from cross-library on-target sgRNAs aggregation approach. This library comprises of 657 on-target sgRNAs targeting 176 consistently underrepresented genes.

**`supp_table_1d_520_genes.xlsx`**
The comprehensive information of 520 genes consistently underrepresented across five CRISPR-KO libraries after library refinement from multi-targeting and off-targeting sgRNAs.

---

## GuideRefine_output_5_libraries/

This subfolder contains the outputs of the GuideRefine quality-filtering pipeline applied to each of the five libraries. For each library, there are four files:

| File type | What it contains |
|---|---|
| `*_refined.tsv` | The final list of retained  sgRNAs that passed all quality filters |
| `*_disposed_sgRNAs.tsv` | sgRNAs that were removed, grouped by the reason they were excluded |
| `*_full_report.xlsx` | A full per-guide annotation table with all scoring and filtering details |
| `*.html` | An interactive quality-control report for the library (open in a web browser) |

The guide RNAs removed during filtering are classified by reason:

| Reason | Explanation |
|---|---|
| Multi-targeting | sgRNA aligns to multiple locations in the genome |
| Single-mismatch off-target | sgRNA align with single mismatch at any location between sgRNA-DNA interface |
| PAM-distal double mismatch | Two mismatches in the 2 bases far away from PAM-distal site (pos 1 & 2 or 19 & 20) |
| No valid alignment | sgRNA does not align anywhere in the genome |
| Non-exonic alignment | sgRNA does not target a protein-coding exon and is unlikely to disrupt gene function |
