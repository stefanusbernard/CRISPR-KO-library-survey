# Missing Files

Scanned: production scripts only (excludes `scripts/backup_scripts/`)

| Script | Line | Issue | Fix |
|--------|------|-------|-----|
| `scripts/00_CRISPR-KO-library-survey-functions.R` | 16 | `source("D:/GitHub/GuideRefine/GuideRefine_functions.R")` — Windows-only path, fails on Mac | Update to `source("../../GuideRefine/GuideRefine_functions.R")` |
| `scripts/Fig02_AB_Fig03A_library_overview.Rmd` | 174 | `data/biogrid_crispr_ko_data/important_data/` directory does not exist; `write_csv` will fail | Run `dir.create("../data/biogrid_crispr_ko_data/important_data", recursive = TRUE)` before knitting |
| `scripts/Fig03_BE_cds_pam_analysis.Rmd` | 203 | `data/biogrid_crispr_ko_data/important_data/removed_genes_all_library.csv` — file missing | Regenerate by running the "Export the removed genes" chunk in `Fig02_AB_Fig03A_library_overview.Rmd` |
| `scripts/data_prep/classify_sgrna.Rmd` | 57–61 | `data/guiderefine_output/T2T-CHM13/old_result_pam_distal_mismatch_only/` directory does not exist | The disposed sgRNA files now live one level up at `data/guiderefine_output/T2T-CHM13/`; update the paths in the T2T-CHM13 chunk |
