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

avana_library_meyers <- import_sgrna_library("~/CRISPR-KO-GuideRefine/public_crispr_library/avana_library.tsv")

annotation <- "~/CRISPR-KO-GuideRefine/annotation_file/T2T-CHM13v2.0_gene_annot_granges.rds"
genome_type_version <- "T2T-CHM13"
bsgenome <- BSgenome.Hsapiens.NCBI.T2TCHM13v2.0

# T2T-CHM13
avana_alignment_meyers <- import_sgrna_library_alignment("~/CRISPR-KO-GuideRefine/object_intermediate/T2T-CHM13_bsgenome/avana_library_aln.csv", "T2T-CHM13")

# import predicted paralog pairs from ENSEMBL 115
ensembl_paralog <- ensembl_paralog_gene_pairs()

paralog_pairs <- ensembl_paralog[[1]]
paralog_individual <- ensembl_paralog[[2]]


avana_classification <- off_target_classification(avana_alignment_meyers)
avana_stratify_multi_target <- multi_target_alignment_bin(avana_classification[[1]], c("perfect", "multi-target guides", "non-targeting"))
list_guides <- avana_stratify_multi_target[[1]] %>% filter(alignment_bin %in% c(2,3,4,5,"> 5")) %>% pull(sgRNA)
avana_annotated <- annotate_all_sgrna(avana_alignment_meyers,                                       
                                      annotation,
                                      bsgenome,
                                      genome_type_version,
                                      list_guides,
                                      paralog_individual,
                                      paralog_pairs)
dbl_target_avana <- avana_stratify_multi_target[[1]] %>% filter(alignment_bin == 2) %>% pull(sgRNA)
avana_multi_paralog <- find_multi_target_paralog(avana_annotated, dbl_target_avana, "avana")


aln_paralog_clean <- avana_annotated %>%
  filter(sgrna %in% dbl_target_avana) %>%
  # select guides without any additional single or double-mismatches
  group_by(sgrna) %>%
  filter(all(n_mismatches == 0)) %>%
  distinct() %>%
  ungroup() %>%
  arrange(sgrna) %>%
  select(sgrna, target) %>%
  group_by(sgrna) %>%
  summarise(target_combined = paste(target, collapse = "+"), .groups = "drop")

# list guides that align with additional single or double-mismatches (we dont care about this)
list_excluded_guides <- avana_annotated %>%
  filter(sgrna %in% dbl_target_avana) %>%
  # select guides with any additional single or double-mismatches
  group_by(sgrna) %>%
  filter(!all(n_mismatches == 0)) %>%
  pull(sgrna) %>%
  unique()

list_paralog_guides <- aln_paralog_clean %>% 
  filter(target_combined %in% c("on-target paralog+off-target paralog_paralog", 
                                "off-target paralog_paralog+on-target paralog",
                                "off-target paralog_paralog+off-target paralog_paralog")) %>% 
  pull(sgrna) %>% 
  unique()

list_non_coding_guides <- aln_paralog_clean %>%
  filter(str_detect(target_combined, "non-coding")) %>%
  pull(sgrna) %>% 
  unique()

list_rest_of_guides <- setdiff(setdiff(unique(aln_paralog_clean$sgrna), list_paralog_guides), list_non_coding_guides)













# if you want to look for paralog + paralog unrelated (not pairs), kindly enable list_paralog_unrelated_guides put it in the summary_df
list_paralog_unrelated_guides <- avana_multi_paralog[[2]] %>% filter(target == "off-target paralog_paralog_unrelated") %>% pull(sgrna) %>% unique()

list_singleton_guides <- avana_multi_paralog[[2]] %>% filter(target %in% c("off-target singleton_paralog", "off-target paralog_singleton", "off-target singleton_singleton")) %>% pull(sgrna) %>% unique()
list_paralog_other_genes <- c(list_paralog_unrelated_guides, list_singleton_guides)

list_non_coding_guides <- avana_multi_paralog[[2]] %>% filter(target %in% c("off-target singleton_non-coding", "off-target paralog_non-coding")) %>% pull(sgrna) %>% unique()


combined_guides <- unique(c(list_paralog_guides, list_paralog_unrelated_guides, list_paralog_other_genes, list_non_coding_guides))


intersect(list_paralog_guides, list_non_coding_guides)
intersect(list_singleton_guides, list_non_coding_guides)


