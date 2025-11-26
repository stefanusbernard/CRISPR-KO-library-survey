library(tidyverse)
library(ggpubr)
library(BSgenome.Hsapiens.NCBI.T2TCHM13v2.0)

source("./00_CRISPR-KO-library-survey-functions.R")
source('~/CRISPR-KO-GuideRefine/GuideRefine_functions.R')

pick_alignment_type <- function(classification_data, alignment_type) {
  data <- read_csv(classification_data) %>%
    filter(alignment %in% alignment_type)
  
  return(data)
}

avana_classification <- pick_alignment_type("../data/sgrna_lfc_data/avana_data/avana_sgrna_alignment_classification.csv", c("perfect", "multi-target guides", "single mismatch", "pam-distal single mismatch"))




obtain_required_data <- function(library, normalized_lfc_data, stratification_data) {
  
  library <- read_tsv(library, col_names = FALSE)
  colnames(library) <- c("sgRNA", "spacer", "gene")
  
  normalized_lfc_data <- read_csv(normalized_lfc_data)
  library_stratification <- stratification_data %>% select(sgRNA, alignment, num_alignments)
  
  # left join
  library_lfc_data <- normalized_lfc_data %>%
    left_join(library, by = c("sgRNA", "spacer", "gene")) %>%
    relocate(sgRNA, .before = spacer) %>%
    relocate(gene, .after = spacer) %>%
    arrange(sgRNA)
  
  terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase")
  # lacZ and luciferase are the control in Toronto V3 library
  terms_for_filtering <- paste(terms, collapse = "|")
  
  # separate control from the library
  
  control_lfc_data <- library_lfc_data %>% filter(str_detect(gene, terms_for_filtering))
  
  library_lfc <- library_lfc_data %>%
    filter(!str_detect(gene, terms_for_filtering)) %>%
    left_join(library_stratification, join_by(sgRNA))
}

avana_library_and_lfc_data <- obtain_required_data("../data/library_data/avana_library.tsv",
                                                   "../data/sgrna_lfc_data/output_normalized/avana_normalized_lfc.csv", 
                                                   avana_classification)



count <- combined_data %>%
  select(sgRNA, alignment, library) %>%
  mutate(alignment = factor(alignment, levels = c("perfect", "single mismatch", "pam-distal single mismatch", "multi-target guides"))) %>%
  distinct() %>%
  group_by(library) %>%
  count(alignment) %>%
  ungroup()
