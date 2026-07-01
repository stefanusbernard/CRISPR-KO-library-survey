# Backup of stale/unused functions from CRISPR-KO-library-survey-functions.R
# Backed up: 1 July 2026
# These functions were identified as unused and removed from the main functions file.

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Function to obtain the list of paralog pairs from Barbara's dataset
# Superseded by ensembl_paralog_gene_pairs()

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

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Calculate the percentage of multi-target, single-mismatch, and PAM-distal double mismatch sgRNA in each libraries
# (from FUNCTION FOR 01_analysis_off-target-sgrna.rmd)

calculate_percentage <- function(library_report, discarded_sgrna, original_sgrna) {
  percentage_discarded <- (sum(library_report[[discarded_sgrna]])/sum(library_report[[original_sgrna]])) * 100
  return(percentage_discarded)
}

# ----------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION FOR 02_analysis_removed_genes.rmd

update_gene_symbol <- function(alignment_data, genes_list) {

  alignment_data <- alignment_data %>%
    mutate(checkGeneSymbols(gene)) %>%
    dplyr::select(-c(x, Approved)) %>%
    dplyr::rename("gene_hgnc" = "Suggested.Symbol") %>%
    relocate(gene_hgnc, .after = gene) %>%
    filter(gene_hgnc %in% genes_list)

  return(alignment_data)
}

# ----------------------------------------------------------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------------------------------------------------------

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
  total <- sum(length(list_excluded_guides), length(list_paralog_guides),
               length(list_non_coding_guides), length(list_rest_of_guides))
  percentage_excluded      <- round(length(list_excluded_guides)   / total * 100, 2)
  percentage_paralog       <- round(length(list_paralog_guides)    / total * 100, 2)
  percentage_non_coding    <- round(length(list_non_coding_guides) / total * 100, 2)
  percentage_rest_of_guides <- round(length(list_rest_of_guides)   / total * 100, 2)

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
