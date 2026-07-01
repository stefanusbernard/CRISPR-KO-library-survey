library(here)
library(ggbreak)
library(HGNChelper)
library(BSgenome)
library(tidyverse)
library(readxl)
library(grid)

source(here::here("scripts", "GuideRefine-functions.R"))

# FUNCTION to import paralog data

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Function to obtain the list of paralog pairs from ENSEMBL 115

ensembl_paralog_gene_pairs <- function(dataset_ensembl = "../data/paralog_data/ensembl_115_human_paralogs_2_dec.txt", 
                                       dataset_hgnc = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt") {
  
  # retrieve ensembl_paralog from ENSEMBL 115
  
  ensembl_paralog <- read_tsv(dataset_ensembl) %>%
    drop_na() %>%
    mutate(gene_pair = paste(`Gene name`, "_", `Human paralogue associated gene name`, sep = "")) %>%
    dplyr::rename("ensembl_gene_id" = `Gene stable ID`,
                  "paralog_name" = `Human paralogue associated gene name`,
                  "percent_identity_human_gene_identical_to_query" = `Paralogue %id. target Human gene identical to query gene`,
                  "percent_identity_query_gene_identical_to_human" = `Paralogue %id. query gene identical to target Human gene`) %>%
    dplyr::select(-`Gene stable ID version`)
  
  gene_pair <- unique(ensembl_paralog$gene_pair)
  gene_list <- unique(ensembl_paralog$`Gene name`)
  
  # HGNC data, select protein-coding gene that has ensembl_gene_id
  # barbara filtered the data by using locus_type == "gene with protein product", our approach is using locus_group
  # if we select locus_group == "protein-coding gene", the locus_type is unique for "gene with protein product"
  
  hgnc_protein_coding <- read_tsv(dataset_hgnc) %>%
    filter(locus_group == "protein-coding gene" & !is.na(ensembl_gene_id)) %>%
    dplyr::rename("gene" = "symbol") %>%
    dplyr::select(gene, locus_group, ensembl_gene_id)
    
  # ENSEMBL gene ID that has paralogs
  
  ensembl_paralog <- ensembl_paralog %>%
    mutate(has_paralogs = TRUE) %>%
    left_join(hgnc_protein_coding, join_by(ensembl_gene_id)) %>%
    relocate(gene, .before = ensembl_gene_id) %>%
    relocate(locus_group, .after = gene) %>%
    dplyr::select(gene, locus_group, ensembl_gene_id, has_paralogs, paralog_name, percent_identity_human_gene_identical_to_query, percent_identity_query_gene_identical_to_human) %>%
    distinct()
  
  hgnc_with_paralogs <- hgnc_protein_coding %>%
    left_join(ensembl_paralog, join_by(ensembl_gene_id, gene, locus_group)) %>%
    mutate(has_paralogs = replace_na(has_paralogs, FALSE))

  ensembl_singleton <- hgnc_with_paralogs %>%
    filter(has_paralogs == FALSE)

  count_paralog_singleton <- hgnc_with_paralogs %>%
    dplyr::select(gene, has_paralogs) %>%
    distinct() %>%
    count(has_paralogs) %>%
    dplyr::rename("number_of_genes" = "n") %>%
    mutate(percentage = round(number_of_genes / sum(number_of_genes) * 100, 2))
  
  
  return(list(gene_pair, gene_list, ensembl_singleton, ensembl_paralog, count_paralog_singleton))
}

# FUNCTION FOR 01_analysis_off-target-sgrna.rmd

# ----------------------------------------------------------------------------------------------------------------------------------------------

import_and_process_report <- function(file_dir) {

  terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase")
  terms_for_filtering <- paste(terms, collapse = "|")

  report <- read_excel(file_dir)
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
  
  # only pick one from duplicate gene original (not gene_hgnc)
  # the reason is the total sgrna from duplicate gene original seems like the total from duplicated data due to small bug in updating gene symbol
  
  check_duplicate <- report %>%
    filter(duplicated(gene) | duplicated(gene, fromLast = TRUE))
  
  print(paste0("Warning!, there are: ", length(unique(check_duplicate$gene_hgnc)), " genes with duplicate rows"))
  
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
  
  # total percentage of failed genes only counted for protein-coding genes defined by HGNC
  percentage_failed_genes <- round(length(unique(failed_genes)) / (length(unique(failed_genes)) + length(unique(passed_genes))), 4) * 100
  
  # because there is additional sgRNA (corrected), hence the formula is:
  # ((original sgRNA (without symbol change) + additional sgRNA (corrected)) - actual total sgRNA) / original sgRNA number (without symbol change)
  percentage_failed_sgrna <-  round(((sum(report$`sgRNA number`) + sum(report$`Additional sgRNA (Corrected)`)) - sum(report$`actual total sgRNA`)) / sum(report$`sgRNA number`), 4) * 100
  
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





# FUNCTION FOR 03_analysis_paralog_off-target-sgrna.rmd

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

  terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase")
  terms_for_filtering <- paste(terms, collapse = "|")

  if (str_detect(genome_type, "T2T-CHM13")) {
    alignment <- read_csv(data_dir) %>% distinct()
    alignment <- alignment %>% filter(!str_detect(gene, terms_for_filtering))
    return(alignment)

  } else if (str_detect(genome_type, "hg38")) {
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

# FUNCTION TO CLASSIFY ALL OFF-TARGET SGRNAS

off_target_classification <- function(library_alignment, pam_distal_single_mismatch, pam_distal_double_mismatch) {
  
  library_alignment_non_targeting <- library_alignment %>%
    filter(is.na(n_mismatches)) %>%
    dplyr::select(sgRNA, spacer, gene, n_mismatches) %>%
    mutate(num_alignments = 0,
           alignment = "non-targeting")
  
  
  classification_alignment <- library_alignment %>%
    group_by(sgRNA, spacer, gene, n_mismatches, .drop = FALSE) %>%
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

  } else {
    stop(paste("Unsupported alignment_type:", paste(alignment_type, collapse = ", ")))
  }


  alignment_bin_df <- count_alignment_type %>% 
    count(alignment_bin) %>%
    dplyr::rename("num_of_alignments" = "alignment_bin",
                  "number_of_sgrnas" = "n")
  
  return(list(count_alignment_type, alignment_bin_df))
  
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

# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION FOR 04_analysis_lfc_off-target-sgrna.rmd

# Reads library TSV + LFC CSV, joins them, and removes control sgRNAs.
# Returns the raw joined data without any alignment stratification applied.
load_library_lfc_raw <- function(library_path, lfc_path) {
  library_df <- read_tsv(library_path, col_names = FALSE)
  colnames(library_df) <- c("sgRNA", "spacer", "gene")

  normalized_lfc_data <- read_csv(lfc_path)

  terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase")
  terms_for_filtering <- paste(terms, collapse = "|")

  normalized_lfc_data %>%
    left_join(library_df, by = c("sgRNA", "spacer", "gene")) %>%
    relocate(sgRNA, .before = spacer) %>%
    relocate(gene,  .after  = spacer) %>%
    arrange(sgRNA) %>%
    filter(!str_detect(gene, terms_for_filtering))
}

# Joins pre-loaded library/LFC data with alignment stratification output from
# count_alignment_bin()[[1]], then drops rows with missing alignment info.
apply_stratification <- function(library_lfc_raw, stratification_data) {
  library_stratification <- stratification_data %>%
    dplyr::select(sgRNA, alignment, num_alignments, alignment_bin)

  library_lfc_raw %>%
    left_join(library_stratification, join_by(sgRNA)) %>%
    drop_na()
}

# Convenience wrapper kept for backwards compatibility with other scripts.
obtain_required_data <- function(library, normalized_lfc_data, stratification_data) {
  raw <- load_library_lfc_raw(library, normalized_lfc_data)
  apply_stratification(raw, stratification_data)
}


calculate_wilcoxon_cles <- function(input_data, which_guides, which_guides_2) {
  
  input_data$alignment <- factor(input_data$alignment, levels = c(which_guides, which_guides_2))
  res <- cliff.delta(mean ~ alignment, data = input_data)
  delta <- res$estimate
  delta
  # cles <- (delta + 1) / 2
  # cles
}



visualize_two_group <- function(library_lfc_data, boxplot_output_name,
                                comparison_group = "multi-target guides",
                                comparison_label = "Multi-target") {

  facet_labels <- c("brunello"   = "Brunello",
                    "toronto_v3" = "TKOv3",
                    "yusa"       = "Yusa",
                    "avana"      = "Avana",
                    "jacquere"   = "Jacquere")

  label_map <- c("perfect" = "On-target")
  label_map[comparison_group] <- comparison_label

  boxplot_lfc <- ggboxplot(library_lfc_data,
                           x         = "alignment",
                           y         = "mean",
                           fill      = "alignment",
                           facet.by  = "library",
                           outliers  = FALSE) +
    stat_compare_means(comparisons  = list(c("perfect", comparison_group)),
                       method       = "wilcox.test",
                       method.args  = list(alternative = "two.sided"),
                       size         = 5,
                       step.increase = 0.05,
                       tip.length   = 0.01,
                       label        = "p.signif",
                       label.y      = 1) +
    scale_y_continuous(limits = c(-3, 2.5), breaks = seq(-2, 2, by = 1)) +
    scale_fill_manual(
      name   = "sgRNA: ",
      values = c("perfect"                    = "#0072B2",
                 "single mismatch"            = "#56B4E9",
                 "pam-distal single mismatch" = "mistyrose",
                 "pam-distal double mismatch" = "#E7F7D5",
                 "multi-target guides"        = "#E69F00"),
      labels = label_map) +
    labs(title = "",
         x     = "",
         y     = "Averaged median sgRNA Log2FC",
         fill  = "Alignment:") +
    theme(plot.margin      = margin(10, 10, 25, 10),
          strip.text.x     = element_text(size = 14),
          axis.text.x      = element_blank(),
          axis.text.y      = element_text(size = 14, colour = "black"),
          axis.ticks.y     = element_line(size = 1.5),
          axis.ticks.x     = element_blank(),
          axis.ticks.length = unit(0.1, "cm"),
          axis.title.y     = element_text(size = 14, vjust = 0.5, margin = margin(r = 20)),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_blank(),
          axis.line.x      = element_line(size = 0.5),
          axis.line.y      = element_line(size = 0.5),
          line             = element_line((size = 2), colour = "black"),
          legend.position  = "top",
          legend.text      = element_text(size = 12)) +
    facet_wrap(~library, ncol = 5, labeller = labeller(library = facet_labels))

  b     <- ggplot_build(boxplot_lfc)
  stats <- b$data[[1]]

  count <- library_lfc_data %>%
    select(sgRNA, alignment, library) %>%
    mutate(alignment = factor(alignment, levels = c("perfect", comparison_group))) %>%
    distinct() %>%
    group_by(library) %>%
    count(alignment) %>%
    ungroup() %>%
    bind_cols(stats %>% select(x, ymin))

  boxplot_lfc <- boxplot_lfc +
    geom_text(
      data  = count,
      aes(x = x, y = ymin - 0.5, label = paste0("n= ", n)),
      angle = 45,
      vjust = 0,
      hjust = 0.5,
      size  = 4.25
    )

  ggsave(
    width    = 7.5,
    height   = 5,
    dpi      = 1000,
    plot     = boxplot_lfc,
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
      axis.text.x = element_text(size = 17, vjust = 0.7, colour = "black"),
      axis.text.y = element_text(size = 17, colour = "black"),
      axis.title.x = element_text(size = 17, margin = margin(t = 10)),
      axis.title.y = element_text(size = 17, margin = margin(r = 10)),
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
      size = 5,
      angle = 45,
      vjust = -0.2,
      hjust = -0.1
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 20, 10, 10))
  
  ggsave(
    width = 5,
    height = 5,
    dpi = 1000,
    plot = boxplot_lfc,
    filename = boxplot_output_name)
  
}

# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION FOR 05_analysis_brunello_jacquere_hits.rmd

fix_gene_symbols <- function(genes) {
  
  checkGeneSymbols(genes) %>%
    as_tibble() %>%
    
    # Use Suggested.Symbol if present, otherwise fall back to original
    mutate(
      Suggested.Symbol = if_else(is.na(Suggested.Symbol),x,Suggested.Symbol)) %>%
    
    # Split ambiguous symbols like "ABCD /// DEFG" into separate rows
    separate_rows(
      Suggested.Symbol,
      sep = "\\s*///\\s*"
    ) %>%
    
    # Optional: clean whitespace just in case
    mutate(
      Suggested.Symbol = str_trim(Suggested.Symbol)
    )
}






