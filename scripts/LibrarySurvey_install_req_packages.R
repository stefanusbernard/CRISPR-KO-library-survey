# Install all packages required for scripts:
#   00_CRISPR-KO-library-survey-functions.R
#   Fig02_AB_Fig03A_library_overview.Rmd
#   Fig02_CJ_lfc_stratification.Rmd
#   Fig03_BE_cds_pam_analysis.Rmd

# BiocManager is required to install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

list_packages_bioconductor <- c(
  # Genome sequences
  "BSgenome",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0",
  # Annotation and genomic features (Fig03_BE)
  "GenomicFeatures",
  "GenomeInfoDbData",
  "txdbmaker",
  # CRISPR design (Fig03_BE)
  "crisprVerse",
  "crisprBase"
)

for (pkg in list_packages_bioconductor) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = TRUE)
  }
}

list_packages_cran <- c(
  # Core data wrangling and plotting
  "tidyverse",
  "readxl",
  "openxlsx",
  # Visualisation utilities
  "ggvenn",
  "ggbreak",
  "VennDiagram",
  "paletteer",
  "RColorBrewer",
  "UpSetR",
  "ggpubr",
  "ggsignif",
  # Gene symbol and enrichment
  "HGNChelper",
  "gprofiler2",
  # Statistics
  "clinfun",
  "rstatix",
  "effsize",
  # Project-root-aware paths (cross-platform)
  "here",
  # Rmd rendering
  "rmarkdown",
  "knitr",
  "pandoc"
)

new_packages <- list_packages_cran[!(list_packages_cran %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Ensure pandoc (system tool required by rmarkdown) is available
if (!rmarkdown::pandoc_available()) {
  message("Pandoc not found. Installing via the pandoc R package...")
  pandoc::pandoc_install()
  pandoc::pandoc_activate()
}

message("Pandoc version: ", rmarkdown::pandoc_version())
