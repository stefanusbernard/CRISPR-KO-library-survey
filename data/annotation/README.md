# Annotation files

This folder contains three genome annotation files required by the analysis scripts. None of them are git-tracked (too large); download them manually before running the pipeline.

---

## Files

### `T2T-CHM13v2.0_gene_annot_granges.rds`

A pre-built R `GRanges` object of T2T-CHM13v2.0 gene annotations (~1.4 MB, **git-tracked**).

**Used by:** `data_preparation_scripts/sgrna_off_target_classification.Rmd` to annotate sgRNA alignment positions.

This file is already in the repository and does not need to be downloaded.

---

### `CCDS.20221027.txt`

NCBI Consensus CDS (CCDS) annotation for hg38 (~10 MB, **not git-tracked**).

**Download from:** [https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/](https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/)

Select the file `CCDS.20221027.txt`.

**Used by:** the hg38 classification chunk in `data_preparation_scripts/sgrna_off_target_classification.Rmd`.

---

### `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz`

Full genome annotation GFF for T2T-CHM13v2.0 (~76 MB compressed, **not git-tracked**).

**Download from:** NCBI — search for assembly `GCF_009914755.1` at [https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) and download the GFF annotation file.

Alternatively via FTP:
```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
```

**Used by:** the `create-txdb` chunk in `scripts/Fig03_BE_cds_pam_analysis.Rmd` (marked `eval=FALSE`). This chunk only needs to be run **once** to build the TxDb SQLite database saved at `data/T2T_data/TxDb_T2T-CHM13v2.sqlite`. Once that database exists, the GFF file is no longer needed.

To build the TxDb manually:

```r
library(txdbmaker)

txdb <- makeTxDbFromGFF(
  file = "data/annotation/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
  format = "gff3"
)

saveDb(txdb, file = "data/T2T_data/TxDb_T2T-CHM13v2.sqlite")
```

---

## Summary

| File | Size | Git-tracked | Download needed |
|---|---|---|---|
| `T2T-CHM13v2.0_gene_annot_granges.rds` | ~1.4 MB | yes | no |
| `CCDS.20221027.txt` | ~10 MB | no | yes (hg38 analysis only) |
| `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz` | ~76 MB | no | yes (one-time TxDb build) |
