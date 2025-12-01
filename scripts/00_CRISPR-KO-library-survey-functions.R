library(ggvenn)
library(ggbreak)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(tidyverse)
library(readxl)
library(VennDiagram)
library(grid)
library(paletteer)
library(RColorBrewer)


source('~/CRISPR-KO-GuideRefine/GuideRefine_functions.R')

# FUNCTION to import paralog data

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Function to obtain the list of paralog pairs from Barbara's dataset

barbara_paralog_gene_pairs <- function(dataset = "../data/paralog_data/barbara_36K_paralog_pairs.csv") {
  
  predicted_paralog_dataset <- read_csv(dataset) %>%
    mutate(sorted_gene_pair_2 = paste(A2, "_", A1, sep = "")) %>%
    relocate(.after = sorted_gene_pair)
  
  gene_pair_1 <- predicted_paralog_dataset$sorted_gene_pair
  gene_pair_2 <- predicted_paralog_dataset$sorted_gene_pair_2
  
  combined <- unique(union(predicted_paralog_dataset$A1, predicted_paralog_dataset$A2))
  
  gene_pair <- unique(union(gene_pair_1, gene_pair_2))
  
  return(list(gene_pair, combined))
}


# Function to obtain the list of paralog pairs from ENSEMBL 115

ensembl_paralog_gene_pairs <- function(dataset_ensembl = "../data/paralog_data/ensembl_115_human_paralogs.txt", 
                                       dataset_hgnc = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt") {
  
  # retrieve ensembl_paralog from ENSEMBL 115
  
  ensembl_paralog <- read_csv(dataset_ensembl) %>%
    drop_na() %>%
    mutate(gene_pair = paste(`Gene name`, "_", `Human paralogue associated gene name`, sep = "")) %>%
    dplyr::rename(ensembl_gene_id = `Gene stable ID`) %>%
    dplyr::select(-`Gene stable ID version`)
  
  gene_pair <- unique(ensembl_paralog$gene_pair)
  gene_list <- unique(ensembl_paralog$`Gene name`)
  
  # HGNC data, select protein-coding gene that has ensembl_gene_id
  # barbara filtered the data by using locus_type == "gene with protein product", our approach is using locus_group
  # if we select locus_group == "protein-coding gene", the locus_type is unique for "gene with protein product"
  
  hgnc <- read_tsv(dataset_hgnc)
  
  hgnc_protein_coding <- hgnc %>%
    filter(locus_group == "protein-coding gene" & !is.na(ensembl_gene_id)) %>%
    dplyr::select(symbol, locus_group, ensembl_gene_id)
  
  # ENSEMBL gene ID that has paralogs
  
  ensembl_paralog <- ensembl_paralog %>%
    mutate(has_paralogs = TRUE) %>%
    dplyr::select(ensembl_gene_id, has_paralogs) %>%
    distinct()
  
  ensembl_singleton <- hgnc_protein_coding %>%
    left_join(ensembl_paralog, join_by(ensembl_gene_id)) %>%
    mutate(has_paralogs = replace_na(has_paralogs, FALSE)) %>%
    filter(has_paralogs == FALSE)
  
  count_paralog_singleton <- hgnc_protein_coding %>%
    left_join(ensembl_paralog, join_by(ensembl_gene_id)) %>%
    mutate(has_paralogs = replace_na(has_paralogs, FALSE)) %>%
    count(has_paralogs) %>%
    dplyr::rename("number_of_genes" = "n") %>%
    mutate(percentage = round(number_of_genes / sum(number_of_genes) * 100, 2))
  
  
  
  return(list(gene_pair, gene_list, ensembl_singleton, ensembl_paralog, count_paralog_singleton))
}

# FUNCTION FOR 01_analysis_off-target-sgrna.rmd

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Calculate the percentage of multi-target, single-mismatch, and PAM-distal double mismatch sgRNA in each libraries

calculate_percentage <- function(library_report, discarded_sgrna, original_sgrna) {
  percentage_discarded <- (sum(library_report[[discarded_sgrna]])/sum(library_report[[original_sgrna]])) * 100
  return(percentage_discarded)
}

import_and_process_report <- function(file_dir, control_terms) {
  
  report <- read_excel(paste(file_dir))
  report <- report %>% filter(!str_detect(gene, terms_for_filtering))
  
  # we noted that our pipeline also removed sgRNAs targeting non protein-coding genes like lncRNA and else, but in this study we focus our investigation to the protein-coding genes
  # protein-coding genes are defined as gene with protein product by HGNC
  
  dataset_hgnc <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
  # to match with paralog investigation later, we obtained protein-coding genes defined by HGNC with established ENSEMBL Gene ID
  hgnc_protein_coding <- read_tsv(dataset_hgnc) %>%
    filter(locus_group == "protein-coding gene" & !is.na(ensembl_gene_id)) %>%
    dplyr::select(symbol, locus_group, ensembl_gene_id)
  
  # first we convert the failed and passed genes in the report from all libraries (except not counted/symbol change)
  # by checkGeneSymbols and update it to a suggested approved symbol
  report <- report %>%
    filter(`quality check` %in% c("fail", "pass"))
  
  # and then we filtered the failed and passed genes if they are protein-coding genes defined by HGNC database
  report_check_gene_symbol <- checkGeneSymbols(report$gene) %>% dplyr::rename("gene" = "x")
  
  # left join the report_check_gene_symbol with the report hence we get the approved gene symbol name
  report <- report %>%
    left_join(report_check_gene_symbol, join_by(gene)) %>%
    dplyr::rename("gene_hgnc" = "Suggested.Symbol") %>%
    mutate(gene_hgnc = if_else(is.na(gene_hgnc), gene, gene_hgnc)) %>%
    relocate(Approved, .after = gene) %>%
    relocate(gene_hgnc, .after = Approved) %>%
    distinct() %>%
    # only select the gene that is defined as protein-coding gene by HGNC data
    filter(gene_hgnc %in% hgnc_protein_coding$symbol)
  
  return(report)
  
}

percentages_and_list_genes <- function(report, list_paralog) {
  
  # pull only the gene list for failed and passed genes
  failed_genes <- report %>%
    filter(`quality check` == "fail") %>%
    pull(gene_hgnc)

  passed_genes <- report %>%
    filter(`quality check` == "pass") %>%
    pull(gene_hgnc)
  
  # pull the dataframe
  failed_genes_df <- report %>%
    filter(`quality check` == "fail")
  
  passed_genes_df <- report %>%
    filter(`quality check` == "pass")
  
  # total percentage of failed genes only counted for protein-coding genes defined by HGNC
  percentage_failed_genes <- round(length(unique(failed_genes)) / (length(unique(failed_genes)) + length(unique(passed_genes))), 4) * 100
  
  # total percentage of failed sgRNAs only counted for protein-coding genes defined by HGNC
  percentage_failed_sgrna <-  round((sum(report$`sgRNA number`) - sum(report$`actual total sgRNA`)) / sum(report$`sgRNA number`), 4) * 100
  
  # percentage of failed genes (paralog and non-paralog protein-coding genes)
  failed_paralog <- failed_genes[unlist(failed_genes) %in% unlist(list_paralog)]
  failed_singleton <- failed_genes[!unlist(failed_genes) %in% unlist(list_paralog)]
  
  percentage_failed_paralog <- round(length(unique(failed_paralog)) / (length(unique(failed_genes))), 4) * 100
  percentage_failed_singleton <- round(length(unique(failed_singleton)) / (length(unique(failed_genes))), 4) * 100
  
  return(list(
    percentage_failed_paralog = percentage_failed_paralog,
    percentage_failed_singleton = percentage_failed_singleton,
    total_percentage_failed_genes = percentage_failed_genes,
    total_percentage_failed_sgrna = percentage_failed_sgrna,
    failed_genes = failed_genes, 
    passed_genes = passed_genes, 
    failed_paralog = failed_paralog,
    failed_singleton = failed_singleton))
}

# FUNCTION FOR 02_analysis_paralog_off-target-sgrna.rmd

# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO IMPORT LIBRARY, ALIGNMENT, OFF-TARGET GUIDES, ETC

import_sgrna_library <- function(data_dir, terms_for_filtering) {
  library <- read_tsv(data_dir, col_names = FALSE)
  colnames(library) <- c('sgRNA', 'spacer', 'gene')
  
  # remove any non-targeting control or control sgRNA
  
  # detect any control/intron control sgRNA
  library <- library %>% filter(!str_detect(gene, terms_for_filtering))
  
  return(library)
}

# import alignment data and keep alignment for normal chromosome only
import_sgrna_library_alignment <- function(data_dir, genome_type) {
  
  if (str_detect("T2T-CHM13", genome_type)) {
    alignment <- read_csv(data_dir) %>% distinct()
    alignment <- alignment %>% filter(!str_detect(gene, terms_for_filtering))
    return(alignment)
    
  } else if (str_detect("hg38", genome_type)) {
    alignment <- read_csv(data_dir) %>% distinct()
    alignment <- alignment %>% filter(!str_detect(gene, terms_for_filtering))
    
    aln_normal_chr <- alignment_normal_chr(alignment)
    return(aln_normal_chr)
    
  }
}

# import list off target guides (multi-target/single-mismatch/double-mismatch-pam-distal)
import_list_off_target_guides <- function(data_dir) {
  off_target_guides <- read_tsv(data_dir)
  
  return(list(
    multi_target_guides = c(na.omit(off_target_guides$multi_target_guides)),
    single_mismatch_guides = c(na.omit(off_target_guides$single_mismatch_guides)),
    pam_distal_single_mismatch_guides = c(na.omit(off_target_guides$pam_distal_single_mismatch_guides)),
    pam_distal_double_mismatch_guides = c(na.omit(off_target_guides$pam_distal_double_mismatch_guides))
  ))
}

# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO STRATIFY MULTI-TARGET SGRNAS

off_target_classification <- function(library_alignment, pam_distal_single_mismatch, pam_distal_double_mismatch) {
  
  library_alignment_non_targeting <- library_alignment %>%
    filter(is.na(n_mismatches)) %>%
    select(sgRNA, n_mismatches) %>%
    mutate(num_alignments = 0,
           alignment = "non-targeting")
  
  
  classification_alignment <- library_alignment %>%
    group_by(sgRNA, n_mismatches, .drop = FALSE) %>%
    summarise(n = n(), .groups = "drop") %>%
    dplyr::rename("num_alignments" = "n") %>%
    filter(n_mismatches != 2) %>%
    mutate(alignment = case_when(
      is.na(n_mismatches) ~ "non-targeting",
      sgRNA %in% unlist(pam_distal_single_mismatch) ~ "pam-distal single mismatch",
      sgRNA %in% unlist(pam_distal_double_mismatch) ~ "pam-distal double mismatch",
      n_mismatches == 0 & num_alignments > 1 ~ "multi-target guides",
      n_mismatches == 0 ~ "perfect",
      n_mismatches == 1 ~ "single mismatch",
      TRUE ~ "other"
    )) %>%
    # make a priority of alignment type
    mutate(alignment = factor(alignment, levels = c("multi-target guides", "pam-distal double mismatch", "pam-distal single mismatch", "single mismatch", "perfect"))) %>%
    
    # sort the data based on the alignment type (factor) and number of alignments (descending, so the highest value is on top of duplicated sgRNA)
    arrange(sgRNA, alignment, desc(num_alignments)) %>%
    group_by(sgRNA) %>%
    distinct(sgRNA, .keep_all = TRUE) %>%
    full_join(library_alignment_non_targeting) %>%
    arrange(sgRNA) %>%
    ungroup()
  
  return(list(classification_alignment, library_alignment_non_targeting))
}


# the reason why we put count_alignment_bin as separate function is because multi-target guides and single mismatch guides
# are different in terms of alignment bin visualization

count_alignment_bin <- function(classification_alignment, alignment_type) {
  
  if("multi-target guides" %in% alignment_type) {
    
    count_alignment_type <- classification_alignment %>%
      # remove any guides with additional single-mismatches and pam-distal double mismatches
      filter(alignment %in% alignment_type)
    
    # bin the number of alignments to see how many guides targeting two until more than 5 different locations in the genome with perfect match
    count_alignment_type$alignment_bin <- cut(count_alignment_type$num_alignments,
                                              breaks = c(-0.5, 0.5:5.5, Inf),  # -0.5 to 5.5 for 0–5, then Inf for >5
                                              labels = c(as.character(0:5), "> 5"),
                                              right = TRUE)
    
  } else if(any(c("single mismatch", "pam-distal single mismatch") %in% alignment_type)) {
    
    # output dataframe for 04_analysis_lfc_off-target-sgrnas.Rmd
    count_alignment_type <- classification_alignment %>%
      # remove any guides with additional multi-target and pam-distal double mismatches
      filter(alignment %in% alignment_type) %>%
      mutate(single_mismatch_alignment = ifelse(alignment %in% c("single mismatch", "pam-distal single mismatch"), num_alignments, 0))
    
    # bin the number of alignments to see how many guides targeting two until more than 5 different locations in the genome with perfect match
    count_alignment_type$alignment_bin <- cut(count_alignment_type$single_mismatch_alignment,
                                              breaks = c(-0.5, 0.5:5.5, Inf),  # -0.5 to 5.5 for 0–5, then Inf for >5
                                              labels = c(as.character(0:5), "> 5"),
                                              right = TRUE)
    
  }
  
  
  alignment_bin_df <- count_alignment_type %>% 
    count(alignment_bin) %>%
    dplyr::rename("num_of_alignments" = "alignment_bin",
                  "number_of_sgrnas" = "n")
  
  return(list(count_alignment_type, alignment_bin_df))
  
}



visualize_off_target_alignment <- function(alignment_bin_df, xlabel, aln_break_1, aln_break_2) {
  
  ggplot(alignment_bin_df, aes(x = num_of_alignments, y = number_of_sgrnas)) +
    geom_col(fill = "#1E88E5") +
    theme_minimal() +
    scale_y_break(c((alignment_bin_df$number_of_sgrnas[alignment_bin_df$num_of_alignments == aln_break_2] + 1500),
                    (alignment_bin_df$number_of_sgrnas[alignment_bin_df$num_of_alignments == aln_break_1] - 1500))) +
    labs(x = xlabel, y = "Number of sgRNAs") +
    geom_text(aes(label = paste("n=", number_of_sgrnas, sep="")),
              nudge_y = 150,
              nudge_x = 0,
              size = 8) +
    theme(axis.text.x = element_text(size = 25, vjust = 0.7, colour = 'black'),
          axis.text.y = element_text(size = 25, colour = "black"),
          axis.title.x = element_text(size = 25, vjust = 0.5, margin = margin(r = 20), colour = "black"),
          axis.title.y = element_text(size = 25, vjust = 0.5, margin = margin(r = 20), colour = "black"),
          axis.ticks.y = element_line(size = 1.5),
          axis.ticks.length = unit(0.3, "cm"),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), 
          axis.text.y.right = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(size = 0.5, colour = "black"),
          axis.line.y = element_line(size = 0.5, colour = "black"),
          legend.position = "none",
          plot.margin = margin(1, 1, 1, 1),
          legend.spacing.x = unit(1, 'cm')) +
    ylim(0, max(alignment_bin_df$number_of_sgrnas) + 1000)
}


# ----------------------------------------------------------------------------------------------------------------------------------------------

# Function to annotate all sgRNAs (in the end only sgRNA targeting coding genes are obtained; this includes on-target and off-target guides)

annotate_sgrna_coding_genes <- function(alignment_data, annotation_data, bsgenome, bsgenome_text) {
  
  if (grepl("\\.txt$", annotation_data)) {
    
    print(paste("Processing:", annotation_data))
    
    # import ccds data
    ccds <- read_ccds_data(ccds_filename = annotation_data)
    
    # transform gene annotation data
    transformed_gene_annotation <- transform_gene_annotation_ccds(ccds, chosen_genome = bsgenome_text)
    
    # take ccds_exon data from transformed gene annotation
    ccds_exon <- transformed_gene_annotation$ccds_exon
    
    # take granges data from transformed gene annotation
    gene_annot_granges <- transformed_gene_annotation$gene_annot_granges_df
    
    
  } else if (grepl("\\.rds$", annotation_data)) {
    
    print(paste("Processing:", annotation_data))
    
    # import .rds data
    gene_annot_granges <- readRDS(file = annotation_data)
    
  }
  
  modified_alignment <- add_cut_pos_pam_pos(alignment_data)
  
  # NCBI if using T2T-CHM13; UCSC if using hg38
  sgrna_target_granges <- make_granges_from_alignment_data(modified_alignment, bsgenome)
  
  gene_df_overlap <- find_overlaps_gene_annotation_and_alignment(sgrna_target_granges, gene_annot_granges)
  
  return(gene_df_overlap)
}


# Function to annotate all sgRNA (coding and non-coding)
# left join the annotated data (gene_df) with not annotate data (alignment_data)
# to differentiate guides targeting coding and non-coding (other genomic loci)

annotate_all_sgrna <- function(alignment_data, annotation_data, bsgenome, bsgenome_text, list_sgrna, list_paralog_individual, list_paralog_pairs){
  
  alignment_data <- alignment_data %>% filter(sgRNA %in% list_sgrna)
  
  # ----------------- Main code to process the data --------------------------
  
  # annotate the alignment_data (the sgRNA will be annotated based on the protein-coding genes location; sgRNA targeting other genomic loci is not annotated)
  gene_df <- annotate_sgrna_coding_genes(alignment_data, annotation_data, bsgenome, bsgenome_text)
  gene_df <- gene_df %>%
    dplyr::rename('actual_chr' = 'chr',
                  'actual_gene' = 'gene')
  
  # this is the standard alignment data (composed of all sgRNA from the alignment; has not been annotated)
  alignment_data <- add_cut_pos_pam_pos(alignment_data)
  alignment_data <- alignment_data %>%
    select(unique_aln_id, sgRNA, spacer, protospacer, gene, chr, n_mismatches, cut_pos) %>%
    dplyr::rename('sgrna' = 'sgRNA')
  
  # alignment_data left joined with the gene_df; so we know which sgRNA are targeting other genomic loci and which sgRNA target coding genes
  alignment_data_gene_df <- alignment_data %>%
    left_join(gene_df, by = join_by(unique_aln_id, sgrna, spacer, protospacer, cut_pos)) %>%
    select(unique_aln_id, sgrna, spacer, protospacer, gene, chr, n_mismatches, cut_pos, actual_chr, actual_gene) %>%
    relocate(gene, .before = actual_gene) %>%
    relocate(chr, .before = actual_chr) %>%
    mutate(actual_gene = replace_na(actual_gene, "non-coding"),
           target = case_when(gene == actual_gene & !gene %in% list_paralog_individual ~ "on-target singleton",
                              gene != actual_gene & !gene %in% list_paralog_individual & !actual_gene %in% list_paralog_individual & !actual_gene == "non-coding" ~ "off-target singleton_singleton",
                              gene != actual_gene & !gene %in% list_paralog_individual & actual_gene %in% list_paralog_individual ~ "off-target singleton_paralog",
                              gene != actual_gene & !gene %in% list_paralog_individual & actual_gene == "non-coding" ~ "off-target singleton_non-coding",
                              gene == actual_gene & gene %in% list_paralog_individual ~ "on-target paralog",
                              gene != actual_gene & !paste(gene, "_", actual_gene, sep = "") %in% list_paralog_pairs & !actual_gene == "non-coding" ~ "off-target paralog_singleton",
                              
                              # obtain paralog - paralog (non-related)
                              gene != actual_gene & !paste(gene, "_", actual_gene, sep = "") %in% list_paralog_pairs & actual_gene %in% list_paralog_individual & !actual_gene == "non-coding" ~ "off-target paralog_paralog_unrelated",
                              # obtain paralog - paralog (related)
                              gene != actual_gene & paste(gene, "_", actual_gene, sep = "") %in% list_paralog_pairs ~ "off-target paralog_paralog",
                              # obtain paralog - non-coding
                              gene != actual_gene & gene %in% list_paralog_individual & actual_gene == "non-coding" ~ "off-target paralog_non-coding",
                              TRUE ~ "Uncategorized"
           )) %>%
    distinct()
  
  if ("Uncategorized" %in% unique(alignment_data_gene_df$target)) {
    print("Uncategorized is in alignment_data_gene_df; kindly check the classification of sgRNA target location")
  } else {
    print("Uncategorized is not in alignment_data_gene_df")
  }
  
  return(alignment_data_gene_df)
}
  

summarize_paralog_targets <- function(sgRNA, paralogs, library) {
  # extract gene names
  genes <- str_extract(sgRNA, "(?<=sg_?)[A-Za-z0-9]+")
  
  # classify paralog vs non-paralog
  paralog_indicator <- ifelse(genes %in% paralogs, "paralog", "singleton")
  
  df_dbl_target_paralog <- data.frame(sgRNA, genes, paralog_indicator, library)
  
  # summary counts + percentages
  df <- data.frame(category = paralog_indicator) %>%
    count(category) %>%
    dplyr::rename("count" = "n") %>%
    mutate(percentage = round(count / sum(count) * 100, 2),
           library = library)
  
  return(list(df, df_dbl_target_paralog))
}

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Function to identify multi-target alignment targeting paralog pairs, strictly only for 2 protein-coding genes (paralog A1 and paralog A2)

find_multi_target_paralog <- function(gene_df, list_query_guides, library_name) {
  
  # select guides without any additional single or double-mismatches and should target paralog pairs
  aln_paralog_clean <- gene_df %>%
    filter(sgrna %in% list_query_guides) %>%
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
  list_excluded_guides <- gene_df %>%
    filter(sgrna %in% list_query_guides) %>%
    # select guides with any additional single or double-mismatches
    group_by(sgrna) %>%
    filter(!all(n_mismatches == 0)) %>%
    pull(sgrna) %>%
    unique()
  
  # find guides that target 1 is paralog and target 2 is paralog (paralog pairs based on ENSEMBL)
  list_paralog_guides <- aln_paralog_clean %>% 
    filter(target_combined %in% c("on-target paralog+off-target paralog_paralog", 
                                  "off-target paralog_paralog+on-target paralog",
                                  "off-target paralog_paralog+off-target paralog_paralog")) %>% 
    pull(sgrna) %>% 
    unique()
  
  # find all guides targeting non-coding (does not matter paralog or not)
  list_non_coding_guides <- aln_paralog_clean %>%
    filter(str_detect(target_combined, "non-coding")) %>%
    pull(sgrna) %>% 
    unique()
  
  # find the rest of guides (does not matter targeting singleton/paralog, or guides targeting the actual same gene)
  list_rest_of_guides <- setdiff(setdiff(unique(aln_paralog_clean$sgrna), list_paralog_guides), list_non_coding_guides)
  
  # calculate the percentage
  percentage_excluded <- round((length(list_excluded_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_non_coding_guides), length(list_rest_of_guides)) * 100), 2)
  
  percentage_paralog <- round((length(list_paralog_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_non_coding_guides), length(list_rest_of_guides)) * 100), 2)

  percentage_non_coding <- round((length(list_non_coding_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_non_coding_guides), length(list_rest_of_guides)) * 100), 2)
  
  percentage_rest_of_guides <- round((length(list_rest_of_guides) / sum(length(list_excluded_guides), length(list_paralog_guides), length(list_non_coding_guides), length(list_rest_of_guides)) * 100), 2)
  
  summary_df <- data.frame(library_name = library_name,
                           total_dbl_target_sgrna = length(list_query_guides),
                           count_sgrna_excluded = length(list_excluded_guides),
                           count_sgrna_paralog = length(list_paralog_guides), 
                           count_sgrna_non_coding = length(list_non_coding_guides),
                           count_sgrna_other_genes = length(list_rest_of_guides), 
                           
                           percentage_sgrna_excluded = percentage_excluded,
                           percentage_sgrna_paralog = percentage_paralog, 
                           percentage_sgrna_non_coding = percentage_non_coding,
                           percentage_sgrna_other_genes = percentage_rest_of_guides)
  
  return(summary_df)
}

# ----------------------------------------------------------------------------------------------------------------------------------------------

# A function to obtain and process sgRNA library Log-fold-change data 

# TODO: compare this code with the code inside the backup script

obtain_required_data <- function(library, normalized_lfc_data, stratification_data) {
  
  library <- read_tsv(library, col_names = FALSE)
  colnames(library) <- c("sgRNA", "spacer", "gene")
  
  normalized_lfc_data <- read_csv(normalized_lfc_data)
  
  library_stratification <- stratification_data %>% select(sgRNA, alignment, num_alignments, alignment_bin)
  
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
    left_join(library_stratification, join_by(sgRNA)) %>%
    drop_na()
  
  return(library_lfc)
}


calculate_wilcoxon_cles <- function(input_data, which_guides, which_guides_2) {
  
  input_data[[1]]$alignment <- factor(input_data[[1]]$alignment, levels = c(which_guides, which_guides_2))
  wilcox.test(input_data[[1]]$mean ~ input_data[[1]]$alignment, alternative = "less")
  res <- cliff.delta(mean ~ alignment, data = input_data[[1]])
  
  delta <- res$estimate
  cles <- (delta + 1) / 2
  cles
}

# FIXME: fix the statistics not showing up

visualize_two_group <- function(library_lfc_data, boxplot_output_name) {

  facet_labels = c("brunello" = "Brunello",
                   "toronto_v3" = "TKOv3",
                   "yusa" = "Yusa (Project Score)",
                   "avana" = "Avana",
                   "jacquere" = "Jacquere")
  
  boxplot_lfc <- ggboxplot(library_lfc_data, 
                           x = "alignment", 
                           y = "mean",
                           fill = "alignment",
                           facet.by = "library",
                           outliers = FALSE) +
    stat_compare_means(comparisons = list(c("perfect", "single mismatch"), 
                                          c("perfect", "pam-distal single mismatch"), 
                                          c("perfect", "multi-target guides"),
                                          c("single mismatch", "pam-distal single mismatch"),
                                          c("single mismatch", "multi-target guides"),
                                          c("pam-distal single mismatch", "multi-target guides")),
                       method = "wilcox.test",
                       method.args = list(alternative = "two.sided"),
                       size = 2.5,
                       step.increase = 0.05,
                       tip.length = 0.01,
                       label = "p.format",
                       label.y = 1) +
    scale_y_continuous(limits = c(-3, 3), breaks = seq(-2, 2, by = 1)) +
    scale_fill_manual(
      name = "sgRNA: ",
      values = c("perfect" = "#0072B2", 
                 "single mismatch" = "powderblue",
                 "pam-distal single mismatch" = "mistyrose",
                 "multi-target guides" = "#E69F00"),
      labels = c("On-target", "single mismatch", "pam-distal single mismatch", "multi-target guides")) +
    labs(title = "", 
         x = "", 
         y = "Averaged median sgRNA Log2FC",
         fill = "Alignment:") +
    theme(plot.margin = margin(10, 10, 25, 10),
          strip.text.x = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.ticks.y = element_line(size = 1.5),
          axis.ticks.x = element_blank(), 
          axis.ticks.length = unit(0.1, "cm"),
          axis.title.y = element_text(size = 14, vjust = 0.5, margin = margin(r = 20)),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(size = 0.5),
          axis.line.y = element_line(size = 0.5), 
          line = element_line((size = 2), colour = 'black'),
          legend.text = element_text(size = 12)) +
    facet_wrap(~library, ncol = 5, labeller = labeller(library = facet_labels))
  
  # get stats
  
  b <- ggplot_build(boxplot_lfc)
  stats <- b$data[[1]]
  
  label_df <- stats %>%
    transmute(
      x = x,
      y = ymax + 0.2,          # add padding
      label = round(ymax, 2)   # or any custom label
    )
  
  count <- combined_data %>%
    select(sgRNA, alignment, library) %>%
    mutate(alignment = factor(alignment, levels = c("perfect", "single mismatch", "pam-distal single mismatch", "multi-target guides"))) %>%
    distinct() %>%
    group_by(library) %>%
    count(alignment) %>%
    ungroup()
  
  count <- count %>% bind_cols(stats %>% select(x, ymin))
  
  boxplot_lfc <-  boxplot_lfc +
                                geom_text(
                                  data = count,
                                  aes(x = x, y = ymin - 0.5, label = paste0("n= ", n)),
                                  angle = 45,
                                  vjust = 0,
                                  hjust = 0.5,
                                  size = 2.5
                                )
                              
  ggsave(
    path = './',
    width = 10,
    height = 5,
    dpi = 1000,
    plot = boxplot_lfc,
    filename = boxplot_output_name)
  
}

# 2nd boxplot: visualization of strafication (increasing number of alignment tend to reduce cell fitness)
visualize_stratify_alignment <- function(library_lfc_data, label_x, boxplot_output_name) {
  
  annotation_df <- library_lfc_data %>%
    count(alignment_bin)
  
  boxplot_lfc <- library_lfc_data %>%
    ggplot(aes(alignment_bin, mean, group = alignment_bin)) +
    geom_boxplot(outlier.shape = NA, size = 0.5) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 0.5) +
    scale_y_continuous(limits = c(-3, 2), breaks = seq(-3, 2, by = 1)) +
    labs(title = "", 
         x = label_x, 
         y = "Averaged median sgRNA Log2FC") +
    theme(
      axis.text.x = element_text(size = 12, vjust = 0.7, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.ticks = element_line(size = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(size = 0.5),
      axis.line.y = element_line(size = 0.5), 
      line = element_line((size = 2), colour = 'black'))
  
  # find the upper whisker of each boxplot for labelling
  b <- ggplot_build(boxplot_lfc)
  boxplot_whisker <- b$data[[1]]$ymax
  
  boxplot_lfc <- boxplot_lfc +
    geom_text(
      data = annotation_df,
      aes(x = alignment_bin, y = boxplot_whisker + 0.1, label = paste0("n= ", n)),
      size = 3.5,
      angle = 45,
      vjust = -0.2,
      hjust = -0.2
    ) +
    coord_cartesian(clip = "off")
  
  ggsave(
    path = './',
    width = 5,
    height = 5,
    dpi = 1000,
    plot = boxplot_lfc,
    filename = boxplot_output_name)
  
}



  
  
  
  
  
  
  
  







