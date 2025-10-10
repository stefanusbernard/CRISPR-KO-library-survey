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

ensembl_paralog <- ensembl_paralog_gene_pairs()

paralog_pairs <- ensembl_paralog[[1]]
paralog_individual <- ensembl_paralog[[2]]

brunello_alignment <- import_sgrna_library_alignment("~/CRISPR-KO-GuideRefine/object_intermediate/T2T-CHM13_bsgenome/broadgpp-brunello-library-contents_aln.csv", "T2T-CHM13")

gene_df <- annotate_sgrna_coding_genes(brunello_alignment, 
                                       "~/CRISPR-KO-GuideRefine/annotation_file/T2T-CHM13v2.0_gene_annot_granges.rds", 
                                       BSgenome.Hsapiens.NCBI.T2TCHM13v2.0, 
                                       "T2T-CHM13")
gene_df <- gene_df %>%
  dplyr::rename('actual_chr' = 'chr',
                'actual_gene' = 'gene')


alignment_data <- add_cut_pos_pam_pos(brunello_alignment)
alignment_data <- alignment_data %>%
  select(unique_aln_id, sgRNA, spacer, protospacer, gene, chr, n_mismatches, cut_pos) %>%
  dplyr::rename('sgrna' = 'sgRNA')


alignment_data_gene_df <- alignment_data %>%
  left_join(gene_df, by = join_by(unique_aln_id, sgrna, spacer, protospacer, cut_pos)) %>%
  select(unique_aln_id, sgrna, spacer, protospacer, gene, chr, n_mismatches, cut_pos, actual_chr, actual_gene) %>%
  relocate(gene, .before = actual_gene) %>%
  relocate(chr, .before = actual_chr) %>%
  mutate(actual_gene = replace_na(actual_gene, "non-coding"),
         target = case_when(gene == actual_gene & !gene %in% paralog_individual ~ "on-target singleton",
                            gene != actual_gene & !gene %in% paralog_individual & !actual_gene %in% paralog_individual & !actual_gene == "non-coding" ~ "off-target singleton_singleton",
                            gene != actual_gene & !gene %in% paralog_individual & actual_gene %in% paralog_individual ~ "off-target singleton_paralog",
                            gene != actual_gene & !gene %in% paralog_individual & actual_gene == "non-coding" ~ "off-target singleton_non-coding",
                            gene == actual_gene & gene %in% paralog_individual ~ "on-target paralog",
                            gene != actual_gene & !paste(gene, "_", actual_gene, sep = "") %in% paralog_pairs & !actual_gene == "non-coding" ~ "off-target paralog_singleton",
                            gene != actual_gene & paste(gene, "_", actual_gene, sep = "") %in% paralog_pairs ~ "off-target paralog_paralog",
                            gene != actual_gene & gene %in% paralog_individual & actual_gene == "non-coding" ~ "off-target paralog_non-coding",
                            TRUE ~ "Uncategorized"
                            )) %>%
  distinct()




brunello_classification <- off_target_classification(brunello_alignment)
brunello_stratify_multi_target <- off_target_alignment_bin(brunello_classification[[1]], c("perfect", "multi-target guides", "non-targeting"))
dbl_target_brunello <- brunello_stratify_multi_target[[1]] %>% filter(alignment_bin == 2) %>% pull(sgRNA)
brunello_dbl_target_genes_paralog <- summarize_paralog_targets(dbl_target_brunello, paralog_individual, "brunello")[[1]]
list_query_guides_paralog <- summarize_paralog_targets(dbl_target_brunello, paralog_individual, "brunello")[[2]] %>% filter(paralog_indicator == "paralog") %>% pull(sgRNA)


# select guides without any additional single or double-mismatches and should target paralog pairs
aln_paralog_clean <- alignment_data_gene_df %>%
  filter(sgrna %in% list_query_guides_paralog) %>%
  # select guides without any additional single or double-mismatches
  group_by(sgrna) %>%
  filter(all(n_mismatches == 0)) %>%
  distinct()

# list guides that align with additional single or double-mismatches
list_excluded_guides <- alignment_data_gene_df %>%
  filter(sgrna %in% list_query_guides_paralog) %>%
  # select guides with any additional single or double-mismatches
  group_by(sgrna) %>%
  filter(!all(n_mismatches == 0)) %>%
  pull(sgrna) %>%
  unique()
  
# length(unique(aln_multi_target_paralog$sgrna))

# of the 2009 double-target guides that originally target gene with paralog, we selected double-target guides without additional single or double mismatches
# we obtained 654 guides fulfill our criteria. Which 

list_paralog_guides <- aln_paralog_clean %>% filter(target == "off-target paralog_paralog") %>% pull(sgrna) %>% unique()
list_singleton_guides <- aln_paralog_clean %>% filter(target %in% c("off-target singleton_paralog", "off-target paralog_singleton")) %>% pull(sgrna) %>% unique()
list_non_coding_guides <- aln_paralog_clean %>% filter(target %in% c("off-target singleton_non-coding", "off-target paralog_non-coding")) %>% pull(sgrna) %>% unique()

percentage_excluded <- round((length(list_excluded_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_singleton_guides), length(list_non_coding_guides)) * 100), 2)
percentage_paralog <- round((length(list_paralog_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_singleton_guides), length(list_non_coding_guides)) * 100), 2)
percentage_singleton <- round((length(list_singleton_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_singleton_guides), length(list_non_coding_guides)) * 100), 2)
percentage_non_coding <- round((length(list_non_coding_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_singleton_guides), length(list_non_coding_guides)) * 100), 2)

list_combined <- c(list_excluded_guides, list_paralog_guides, list_singleton_guides, list_non_coding_guides)


investigate_guides <- setdiff(list_query_guides_paralog, list_combined)

# seems like these guides are multi-target guides targeting the same gene
df <- alignment_data_gene_df %>%
  filter(sgrna %in% investigate_guides)



count_single_mismatch <- avana_classification[[1]] %>%
  mutate(single_mismatch_alignment = ifelse(alignment == "single mismatch", num_alignments, 0))

# bin the number of alignments to see how many guides targeting two until more than 5 different locations in the genome with perfect match

count_single_mismatch$alignment_bin <- cut(count_single_mismatch$single_mismatch_alignment,
                                                breaks = c(-0.5, 0.5:5.5, Inf),  # -0.5 to 5.5 for 0–5, then Inf for >5
                                                labels = c(as.character(0:5), "> 5"),
                                                right = TRUE)
alignment_bin_df <- count_single_mismatch %>% 
  count(alignment_bin) %>%
  dplyr::rename("num_of_alignments" = "alignment_bin",
                "number_of_sgrnas" = "n")















