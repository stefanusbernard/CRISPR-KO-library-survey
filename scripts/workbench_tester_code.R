library(tidyverse)
library(ggpubr)
library(BSgenome.Hsapiens.NCBI.T2TCHM13v2.0)

source("./00_CRISPR-KO-library-survey-functions.R")
source('~/CRISPR-KO-GuideRefine/GuideRefine_functions.R')


dataset_ensembl <- "../data/paralog_data/ensembl_115_human_paralogs.txt"
dataset_hgnc <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
  
# retrieve ensembl_paralog from ENSEMBL 115

ensembl_paralog <- read_csv(dataset_ensembl) %>%
  drop_na() %>%
  mutate(gene_pair = paste(`Gene name`, "_", `Human paralogue associated gene name`, sep = "")) %>%
  dplyr::rename(ensembl_gene_id = `Gene stable ID`) %>%
  dplyr::select(-`Gene stable ID version`) %>%
  distinct()

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
  
  
