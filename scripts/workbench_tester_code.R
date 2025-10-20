library(tidyverse)
library(VennDiagram)
library(ggvenn)
library(ggbreak)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.T2TCHM13v2.0)

source("./00_CRISPR-KO-library-survey-functions.R")
source('~/CRISPR-KO-GuideRefine/GuideRefine_functions.R')


terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase")
# lacZ and luciferase are the control in Toronto V3 library

terms_for_filtering <- paste(terms, collapse = "|")

disposed_avana <- import_list_off_target_guides("~/CRISPR-KO-GuideRefine/output_cleaning/T2T-CHM13/avana_library_disposed_sgRNAs.tsv")

annotation <- "~/CRISPR-KO-GuideRefine/annotation_file/T2T-CHM13v2.0_gene_annot_granges.rds"
genome_type_version <- "T2T-CHM13"
bsgenome <- BSgenome.Hsapiens.NCBI.T2TCHM13v2.0

# T2T-CHM13
avana_alignment_meyers <- import_sgrna_library_alignment("~/CRISPR-KO-GuideRefine/object_intermediate/T2T-CHM13_bsgenome/avana_library_aln.csv", "T2T-CHM13")
avana_classification <- off_target_classification(avana_alignment_meyers, disposed_avana["pam_distal_double_mismatch_guides"])


avana_single_mismatch <- avana_classification[[1]] %>% filter(alignment == "single mismatch") %>% pull(sgRNA)
avana_single_mismatch_df <- avana_alignment_meyers %>%
  filter(sgRNA %in% avana_single_mismatch & n_mismatches == 1) %>%
  arrange(sgRNA)




