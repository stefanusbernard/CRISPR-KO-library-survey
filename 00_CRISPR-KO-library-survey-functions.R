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


source('../CRISPR-KO-library-refinement-pipeline/sgRNA_library_refinement_pipeline_functions.R')


# FUNCTION FOR 01_Percentage_off-target-sgRNAs.rmd

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Function to obtain the list of paralog pairs from Barbara's dataset

list_paralog_gene_pairs <- function(dataset = "./barbara_36K_paralog_pairs.csv") {
  
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

ensembl_paralog_gene_pairs <- function(dataset = "./ensembl_115_human_paralogs.txt") {
  
  ensembl_paralog <- read_csv(dataset) %>% 
    drop_na() %>%
    mutate(gene_pair = paste(`Gene name`, "_", `Human paralogue associated gene name`, sep = ""))
  
  gene_pair <- unique(ensembl_paralog$gene_pair)
  gene_list <- unique(ensembl_paralog$`Gene name`)
  
  return(list(gene_pair, gene_list))
  
}



# Calculate the percentage of multi-target, single-mismatch, and PAM-distal double mismatch sgRNA in each libraries

calculate_percentage <- function(library_report, discarded_sgrna, original_sgrna) {
  percentage_discarded <- (sum(library_report[[discarded_sgrna]])/sum(library_report[[original_sgrna]])) * 100
  return(percentage_discarded)
}

import_library <- function(lib_dir){
  library <- read_tsv(lib_dir, col_names = FALSE)
  colnames(library) <- c("sgrna", "spacer", "gene")
  
  return(library)
}

percentages_and_list_genes <- function(report, list_paralog) {
  
  # pull only the gene list for failed and passed genes
  failed_genes <- report %>%
    filter(`quality check` == "fail") %>%
    pull(gene)
  
  passed_genes <- report %>%
    filter(`quality check` == "pass") %>%
    pull(gene)
  
  # pull the dataframe
  failed_genes_df <- report %>%
    filter(`quality check` == "fail")
  
  passed_genes_df <- report %>%
    filter(`quality check` == "pass")
  
  # exclude not counted (symbol change) for total percentage of failed sgRNAs
  # cleaned_report <- report %>%
  #   filter(!`quality check` == "not counted (symbol change)")
  
  # total percentage of failed genes
  percentage_failed_genes <- round(length(unique(failed_genes)) / (length(unique(failed_genes)) + length(unique(passed_genes))), 4) * 100
  
  # total percentage of failed sgRNAs (check with this again later)
  percentage_failed_sgrna <-  round((sum(report$`sgRNA number`) - sum(report$`actual total sgRNA`)) / sum(report$`sgRNA number`), 4) * 100
  
  # percentage of failed genes (paralog and non-paralog protein-coding genes)
  failed_paralog <- failed_genes[unlist(failed_genes) %in% unlist(list_paralog)]
  failed_not_paralog <- failed_genes[!unlist(failed_genes) %in% unlist(list_paralog)]
  
  percentage_failed_paralog <- round(length(unique(failed_paralog)) / (length(unique(failed_genes))), 4) * 100
  percentage_failed_not_paralog <- round(length(unique(failed_not_paralog)) / (length(unique(failed_genes))), 4) * 100
  
  return(list(
    percentage_failed_paralog = percentage_failed_paralog,
    percentage_failed_not_paralog = percentage_failed_not_paralog,
    total_percentage_failed_genes = percentage_failed_genes,
    total_percentage_failed_sgrna = percentage_failed_sgrna,
    failed_genes = failed_genes, 
    passed_genes = passed_genes, 
    failed_paralog = failed_paralog,
    failed_not_paralog = failed_not_paralog))
}



# FUNCTION FOR 02_Paralog_off-target-sgRNAs.rmd

# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO IMPORT LIBRARY, ALIGNMENT, OFF-TARGET GUIDES, ETC

import_sgrna_library <- function(data_dir) {
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
    multi_target_guides = na.omit(off_target_guides$multi_target_guides),
    single_mismatch_guides = na.omit(off_target_guides$single_mismatch_guides),
    pam_distal_double_mismatch_guides = na.omit(off_target_guides$pam_distal_double_mismatch_guides)
  ))
}

# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO STRATIFY MULTI-TARGET SGRNAS

off_target_classification <- function(library_alignment) {
  
  # find any non-targeting alignment in the alignment data
  library_alignment_non_targeting <- library_alignment %>%
    filter(is.na(n_mismatches)) %>%
    select(sgRNA, n_mismatches) %>%
    mutate(num_alignments = 0,
           alignment = "non-targeting")
  
  # classify the on-target guides, multi-target guides, and single-mismatch alignment guides
  classification_alignment <- library_alignment %>%
    group_by(sgRNA, n_mismatches, .drop = FALSE) %>%
    summarise(n = n(), .groups = "drop") %>%
    dplyr::rename("num_alignments" = "n") %>%
    filter(n_mismatches != 2) %>%
    mutate(alignment = case_when(
      is.na(n_mismatches) ~ "non-targeting",
      n_mismatches == 0 & num_alignments > 1 ~ 'multi-target guides',
      n_mismatches == 0 ~ 'perfect',
      n_mismatches == 1 ~ 'single mismatch',
      TRUE ~ 'other'
    )) %>%
    group_by(sgRNA) %>%
    mutate(alignment = if (any(alignment == "multi-target guides")) "multi-target guides" else alignment) %>%
    # mutate(alignment = if (any(alignment == "single mismatch")) "single mismatch" else alignment) %>%
    ungroup() %>%
    filter(!(alignment == "multi-target guides" & n_mismatches > 0),
           !(alignment == "single mismatch" & n_mismatches == 0)) %>%
    full_join(library_alignment_non_targeting) %>%
    arrange(sgRNA)
  
  return(list(classification_alignment, library_alignment_non_targeting))
}

off_target_alignment_bin <- function(classification_alignment, alignment_type) {
  
  count_on_target_alignments <- classification_alignment %>%
    filter(alignment %in% alignment_type)
  
  # bin the number of alignments to see how many guides targeting two until more than 8 different locations in the genome with perfect match
  
  count_on_target_alignments$alignment_bin <- cut(count_on_target_alignments$num_alignments,
                                                  breaks = c(-0.5, 0.5:5.5, Inf),  # -0.5 to 5.5 for 0–5, then Inf for >5
                                                  labels = c(as.character(0:5), "> 5"),
                                                  right = TRUE)
  
  alignment_bin_df <- count_on_target_alignments %>% 
    count(alignment_bin) %>%
    dplyr::rename("num_of_alignments" = "alignment_bin",
                  "number_of_sgrnas" = "n")
  
  return(list(count_on_target_alignments, alignment_bin_df))
  
}


visualize_off_target_alignment <- function(alignment_bin_df, xlabel) {
  
  ggplot(alignment_bin_df, aes(x = num_of_alignments, y = number_of_sgrnas, fill = "#E66100")) +
    geom_col() +
    theme_minimal() +
    scale_y_break(c((alignment_bin_df$number_of_sgrnas[alignment_bin_df$num_of_alignments == 2] + 1500),
                    (alignment_bin_df$number_of_sgrnas[alignment_bin_df$num_of_alignments == 1] - 1500))) +
    labs(x = xlabel, y = "Number of sgRNAs") +
    geom_text(aes(label = paste("n=", number_of_sgrnas, sep="")),
              nudge_y = 100,
              nudge_x = 0,
              size = 9) +
    theme(axis.text.x = element_text(size = 25, vjust = 0.7, colour = 'black'),
          axis.text.y = element_text(size = 25, colour = "black"),
          axis.title.x = element_text(size = 25, vjust = 0.5, margin = margin(r = 20), colour = "black"),
          axis.title.y = element_text(size = 25, vjust = 0.5, margin = margin(r = 20), colour = "black"),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(), 
          axis.text.y.right = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(size = 0.5, colour = "black"),
          axis.line.y = element_line(size = 0.5, colour = "black"),
          legend.position = "none",
          plot.margin = margin(1, 1, 1, 1),
          legend.spacing.x = unit(1, 'cm'))
}


# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO ANNOTATE MULTI-TARGET-/OFF-TARGET SGRNAS

# annotate off-target sgRNAs
annotate_off_target <- function(alignment_data, annotation_data, bsgenome, bsgenome_text) {
  
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
  
  # ccds <- read_ccds_data('../CCDS.20221027.txt')
  # transformed_ccds <- transform_gene_annotation_ccds(ccds)
  
  modified_alignment <- add_cut_pos_pam_pos(alignment_data)
  
  # NCBI if using T2T-CHM13; UCSC if using hg38
  off_target_granges <- make_granges_from_alignment_data(modified_alignment, bsgenome)
  
  gene_df_overlap <- find_overlaps_gene_annotation_and_alignment(off_target_granges, gene_annot_granges)
  
  return(gene_df_overlap)
}


# off target annotation (can be multi-target, single-mismatch, or pam-distal double mismatches) used to annotate whether the sgRNA target random genomic loci, protein-coding genes or intended target gene

off_target_annotation <- function(alignment_data, annotation_data, bsgenome, bsgenome_text, library, list_sgrna, num_mismatches, library_name){
  
  alignment_data <- alignment_data %>% filter(sgRNA %in% list_sgrna)
  
  gene_df <- annotate_off_target(alignment_data, annotation_data, bsgenome, bsgenome_text)
  gene_df <- gene_df %>%
    dplyr::rename('actual_chr' = 'chr',
                  'actual_gene' = 'gene')
  
  alignment_data <- add_cut_pos_pam_pos(alignment_data)
  alignment_data <- alignment_data %>%
    select(unique_aln_id, sgRNA, spacer, protospacer, gene, chr, n_mismatches, cut_pos) %>%
    dplyr::rename('sgrna' = 'sgRNA')
  
  # one version of alignment data for paralog analysis (you will need the actual gene to detect which multi-target sgRNA targeting only two protein-coding genes - one is the actual gene, the other one is the paralog pairs)
  alignment_data_gene_df <- alignment_data %>%
    filter(n_mismatches == num_mismatches) %>%
    left_join(gene_df, by = join_by(unique_aln_id, sgrna, spacer, protospacer, cut_pos)) %>%
    select(unique_aln_id, sgrna, spacer, protospacer, gene, chr, n_mismatches, cut_pos, actual_chr, actual_gene) %>%
    mutate(target = case_when(gene == actual_gene ~ 'actual target gene',
                              gene != actual_gene ~ 'other protein-coding genes',
                              is.na(actual_gene) ~ 'other genomic loci')) %>%
    distinct()
  
  summary_target_df <- alignment_data_gene_df %>%
    filter(target != "actual target gene") %>%
    select(sgrna, target) %>%
    group_by(sgrna) %>%
    # first removal of duplicate data to remove duplicate alignment to multiple locations 
    distinct() %>%
    # there may be a case where a single sgRNA can perfectly target both protein-coding gene and genomic loci, this code is used assign those sgRNA those sgRNA as targeting protein-coding gene only (which may result in duplicated data)
    mutate(target = ifelse(n() > 1, 'other protein-coding genes', target)) %>%
    ungroup() %>%
    # second removal of duplicate data
    distinct()
  
  # there is a possibility a single guide target multiple locations of the same gene and other gene/random genomic locations, this code is used to separate the sgRNA that align multiple times to the same gene with those align to other gene/genomic loci
  multiple_aln_same_target <- alignment_data_gene_df %>%
    filter(target == "actual target gene" & !sgrna %in% summary_target_df$sgrna) %>%
    count(sgrna) %>%
    filter(n > 1) %>%
    dplyr::rename("number_of_alignments" = "n")
  
  # count affected genes
  
  sgrna_target_genes <- summary_target_df %>%
    filter(target == "other protein-coding genes") %>%
    distinct() %>%
    pull(sgrna)
  
  affected_genes <- alignment_data_gene_df %>%
    filter(sgrna %in% sgrna_target_genes & target == "other protein-coding genes") %>%
    select(actual_gene) %>%
    distinct() %>%
    pull(actual_gene)
  
  df_affected_genes_percentage <- data.frame(target = c("other genomic loci", "other protein-coding genes"),
                                             count_genes = c(0, length(affected_genes)),
                                             percent_genes = c(0, round(length(affected_genes)/length(unique(library$gene))*100, 2)))
  
  
  # create frequency table
  
  calculate_freq_table <- summary_target_df %>%
    select(target) %>%
    count(target) %>%
    dplyr::rename('count_sgrna' = 'n') %>%
    arrange(desc(count_sgrna)) %>%
    mutate(percent_sgrna_among_library = round(count_sgrna/length(unique(library$sgRNA))*100 ,2),
           percent_sgrna_among_off_target = round(count_sgrna/length(unique(alignment_data_gene_df$sgrna))*100, 2),
           library = library_name) %>%
    relocate(library, .before = target) %>%
    left_join(df_affected_genes_percentage, join_by(target))
  
  
  return(list(
    aln_data_df = alignment_data_gene_df,
    multiple_aln_same_target = multiple_aln_same_target,
    summary_df = summary_target_df, 
    freq_table = calculate_freq_table,
    list_affected_genes = affected_genes))
}


# ----------------------------------------------------------------------------------------------------------------------------------------------


# Function to identify multi-target alignment targeting paralog pairs, strictly only for 2 protein-coding genes (paralog A1 and paralog A2)

find_multi_target_paralog <- function(alignment_gene_df, list_guides, paralog_dataset, library, library_name) {
  dbl_target_classification <- alignment_gene_df %>%
    filter(sgrna %in% list_guides) %>%
    mutate(target = case_when(paste(gene, actual_gene, sep ='_') %in% paralog_dataset ~ 'targeting paralog gene',
                              gene == actual_gene ~ 'actual gene',
                              is.na(actual_gene) ~ 'other genomic loci',
                              TRUE ~ 'other protein-coding gene'))
  
  aln_multi_target_paralog <- dbl_target_classification %>%
    group_by(sgrna) %>%
    filter(all(c("targeting paralog gene", "actual gene") %in% target),
           n_distinct(target) > 1) %>%
    ungroup()
  
  list_paralog_guides <- dbl_target_classification %>% filter(target == "targeting paralog gene") %>% pull(sgrna) %>% unique()
  list_other_genes_guides <- dbl_target_classification %>% filter(target == "other protein-coding gene") %>% pull(sgrna) %>% unique()
  list_other_genomic_loci <- dbl_target_classification %>% filter(target == "other genomic loci") %>% pull(sgrna) %>% unique()
  
  list_other_genomic_loci_guides <- setdiff(
    list_other_genomic_loci,
    c(list_other_genes_guides, list_paralog_guides)
  )
  
  
  
  list_not_paralog_guides <- union(list_other_genes_guides, list_other_genomic_loci_guides)
  
  # obtain the affected genes by looking at the actual target from double-target guides
  list_paralog_genes <- unique(aln_multi_target_paralog$actual_gene)
  
  # count percentage of the double target paralog guides
  
  count_sgrna_paralog <- length(list_paralog_guides)
  count_sgrna_other_genes <- length(list_other_genes_guides)
  count_sgrna_other_genomic_loci <- length(list_other_genomic_loci_guides)
  
  percentage_sgrna_paralog <- round((count_sgrna_paralog/sum(count_sgrna_paralog, count_sgrna_other_genes, count_sgrna_other_genomic_loci) * 100), 2)
  percentage_sgrna_other_genes <- round((count_sgrna_other_genes/sum(count_sgrna_paralog, count_sgrna_other_genes, count_sgrna_other_genomic_loci) * 100), 2)
  percentage_sgrna_other_genomic_loci <- round((count_sgrna_other_genomic_loci/sum(count_sgrna_paralog, count_sgrna_other_genes, count_sgrna_other_genomic_loci) * 100), 2)
  
  summary_df <- data.frame(library_name,
                           count_sgrna_paralog, 
                           count_sgrna_other_genes, 
                           count_sgrna_other_genomic_loci, 
                           percentage_sgrna_paralog, 
                           percentage_sgrna_other_genes, 
                           percentage_sgrna_other_genomic_loci)
  
  # count percentage of genes affected
  
  paralog_genes_affected <- length(list_paralog_genes)
  total_genes_library <- length(unique(library$gene))
  percentage_genes_affected <- round((paralog_genes_affected/total_genes_library)*100, 2)
  
  count_gene <- data.frame(library_name, paralog_genes_affected, total_genes_library, percentage_genes_affected)
  
  return(list(df_target_classification = dbl_target_classification,
              aln_multi_target_paralog_df = aln_multi_target_paralog,
              summary_count_df = summary_df,
              multi_target_paralog_guides = list_paralog_guides,
              multi_target_other_genes = list_other_genes_guides,
              multi_target_other_genomic_loci = list_other_genomic_loci_guides,
              paralog_affected_list = list_paralog_genes,
              paralog_affected_df = count_gene))
  
}







