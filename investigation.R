library(tidyverse)
library(VennDiagram)
library(ggvenn)
library(ggbreak)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.T2TCHM13v2.0)

source('../CRISPR-KO-library-refinement-pipeline/sgRNA_library_refinement_pipeline_functions.R')
source("./00_CRISPR-KO-library-survey-functions.R")




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


predicted_paralog_dataset <- read_csv("./barbara_36K_paralog_pairs.csv") %>%
  mutate(sorted_gene_pair_2 = paste(A2, "_", A1, sep = "")) %>%
  relocate(sorted_gene_pair_2, .after = sorted_gene_pair)
# 
# gene_pair_1 <- predicted_paralog_dataset$sorted_gene_pair
# gene_pair_2 <- predicted_paralog_dataset$sorted_gene_pair_2
# 
# combined <- unique(union(predicted_paralog_dataset$A1, predicted_paralog_dataset$A2))
# 
# gene_pair <- unique(union(gene_pair_1, gene_pair_2))


ensembl_paralog <- read_csv("./ensembl_115_human_paralogs.txt") %>% 
  drop_na() %>%
  mutate(gene_pair = paste(`Gene name`, "_", `Human paralogue associated gene name`, sep = ""))

gene_pair <- unique(ensembl_paralog$gene_pair)
gene_list <- unique(ensembl_paralog$`Gene name`)



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


dbl_target_classification <- avana_multi_protein_coding_genes[[1]] %>%
  # filter(sgrna %in% list_guides) %>%
  mutate(pairs = paste(gene, actual_gene, sep ='_')) %>%
  mutate(target = case_when(paste(gene, actual_gene, sep ='_') %in% gene_pair ~ 'targeting paralog gene',
                            gene == actual_gene ~ 'actual gene',
                            is.na(actual_gene) ~ 'other genomic loci',
                            TRUE ~ 'other protein-coding gene'))


  

  
# TESTING
  
# annotation <- "../CRISPR-KO-library-refinement-pipeline/annotation_file/T2T-CHM13v2.0_gene_annot_granges.rds"
# genome_type_version <- "T2T-CHM13"
# bsgenome <- BSgenome.Hsapiens.NCBI.T2TCHM13v2.0

terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase")
# lacZ and luciferase are the control in Toronto V3 library

terms_for_filtering <- paste(terms, collapse = "|")

# 
# avana_library_meyers <- import_sgrna_library("../CRISPR-KO-library-refinement-pipeline/public_crispr_library/avana_library.tsv")
avana_alignment_meyers <- import_sgrna_library_alignment("../CRISPR-KO-library-refinement-pipeline/object_intermediate/T2T-CHM13_bsgenome/avana_library_aln.csv", "T2T-CHM13")






  
  