# T2T Annotation Scripts

These scripts prepare the T2T-CHM13v2.0 genome resources required by the main analysis pipeline. Run them once before running any figure scripts.

---

## Required input files

Both scripts depend on two files downloaded from the NCBI FTP server. Place them in `data/annotation/` before running.

### 1. T2T-CHM13v2.0 genomic GFF file

Download via FTP or HTTPS:

```bash
# HTTPS (curl)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz

# FTP (wget)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
```

Expected filename: `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz`

### 2. Assembly report (chromosome ↔ RefSeq accession mapping)

```bash
# HTTPS (curl)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt

# FTP (wget)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt
```

Expected filename: `T2T-CHM13v2_assembly_report.txt`

This file maps RefSeq accession numbers (e.g. `NC_060925.1`) to chromosome names (e.g. `chr1`). It is used in `01_T2T_gff_to_ccds_conversion.Rmd` to produce human-readable chromosome labels.

> Both files can also be found by browsing the NCBI assembly page for accession **GCF_009914755.1** at `https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/`.

---

## Scripts

### `01_T2T_gff_to_ccds_conversion.Rmd`

Converts the T2T-CHM13v2.0 GFF annotation into a CCDS-like format and builds the GRanges object used by the sgRNA classification pipeline.

**Requires:** `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz`, `T2T-CHM13v2_assembly_report.txt`

**What it does:**

1. Extracts canonical (RefSeq Select) mRNA transcripts from the GFF using a `zgrep` command-line filter — no need to decompress the full file.
2. Reads all CDS features from the same GFF and joins them to the canonical transcripts.
3. Maps RefSeq accession numbers to chromosome names using the assembly report.
4. Collapses CDS ranges per gene and formats the result into a CCDS-compatible table (`cds_gene`).
5. Unnests per-exon coordinates, assigns exon codes, and builds a `GRanges` object with T2T-CHM13v2.0 seqinfo.
6. Saves the result as `T2T-CHM13v2.0_gene_annot_granges.rds` — the annotation object loaded by `sgrna_off_target_classification.Rmd`.

**Produces:** `data/annotation/T2T-CHM13v2.0_gene_annot_granges.rds`

**R packages required:**

| Package | Purpose |
|---|---|
| tidyverse | Data wrangling |
| rtracklayer | GFF import |
| GenomicRanges | GRanges construction |
| data.table | Fast GFF parsing via `fread` + shell command |
| ape | (utility support) |

---

### `02_forge_T2T_bsgenome.Rmd`

Builds and installs the `BSgenome.Hsapiens.NCBI.T2TCHM13v2.0` package from the NCBI T2T-CHM13v2.0 assembly. Run this **once** on any new machine before running the main analysis scripts.

**Requires:** internet access (downloads ~3 GB genome sequence from NCBI during `forgeBSgenomeDataPkgFromNCBI()`).

**What it does:**

1. Installs `BSgenomeForge` from Bioconductor.
2. Calls `forgeBSgenomeDataPkgFromNCBI()` with assembly accession `GCF_009914755.1` to download the genome sequence and build a local BSgenome source package.
3. Installs the built package from source into the local R library.
4. Loads `BSgenome.Hsapiens.NCBI.T2TCHM13v2.0` to confirm the installation succeeded.

**Produces:** local R package `BSgenome.Hsapiens.NCBI.T2TCHM13v2.0` (installed into the R library)

> The download timeout is set to 2 hours (`options(timeout = 7200)`) because the T2T genome is ~3 GB. Ensure a stable internet connection before running.

**R packages required:**

| Package | Purpose |
|---|---|
| BSgenomeForge | Forge BSgenome packages from NCBI assemblies |
