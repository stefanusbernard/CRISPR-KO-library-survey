library(tidyverse)
library(HGNChelper)
library(openxlsx)
library(MAGeCKFlute)

path = "data/result_gene_hits_change/tkov3_olivieri_paper/"

all_files <- list.files(
  path = path,
  pattern = "_output_.*\\.txt$",
  full.names = TRUE
)

drug_names <- sort(
  unique(
    sub("_output_.*$", "", basename(all_files))
  )
)

pad_to <- function(x, n) {
  length(x) <- n
  x
}

results <- list()
df_list_genes <- list()

for (drug in drug_names) {
  
  print(paste0("analysing ", drug))
  
  original_file = paste0(path, drug, "_output_original.txt")
  refined_file = paste0(path, drug, "_output_refined.txt")
  
  # filter FDR < 0.15 and normZ < -3 or > 6
  
  filter_and_process_data <- function(data) {
    
    data <- read_tsv(data)
    lower_outlier_limit <- -CutoffCalling(data$normZ, scale = 3)
    upper_outlier_limit <- CutoffCalling(data$normZ, scale = 3)
    
    data <- data %>%
      # Olivieri Screen threshold
      # filter((fdr_synth < 0.15 & normZ < -3) | (fdr_supp < 0.15 & normZ > 6)) %>%
      
      # MAGeCK cutoff threshold
      filter((fdr_synth < 0.15 & normZ < lower_outlier_limit) | (fdr_supp < 0.15 & normZ > upper_outlier_limit)) %>%
      pull(GENE)
    
    # pick this if we dont update the gene symbol
    # return(data)
    
    # pick this if we want to update the gene symbol
    # updating the gene symbol causes a differences of the genes hits discovered by Olivieri (they are updated)
    data <- checkGeneSymbols(data)
    return(data$Suggested.Symbol)
    
  }
  
  # all genes
  set_original <- filter_and_process_data(original_file)
  set_refined <- filter_and_process_data(refined_file)
  
  # exclusive genes
  original <- setdiff(set_original, set_refined)
  refined <- setdiff(set_refined, set_original)
  
  # overlap genes
  overlap <- intersect(set_original, set_refined)
  
  max_len <- max(length(overlap), length(original), length(refined))
  
  # make a table list of differences between hits
  results[[length(results) + 1]] <- data.frame(
    drug = drug,
    total_original_hits = length(set_original),
    total_refined_hits = length(set_refined),
    overlap_hits = length(overlap),
    hits_only_in_original = length(original),
    hits_only_in_refined  = length(refined)
  )
  
  # make a table composed a list of genes and output it as excel
  data_list_genes <- data.frame(
    overlap  = pad_to(overlap,  max_len),
    original = pad_to(original,  max_len),
    refined  = pad_to(refined,  max_len),
    stringsAsFactors = FALSE
  )

  names(data_list_genes) <- paste0(
    drug, "_",
    c("overlap_hits", "original_only_hits", "refined_only_hits")
  )
  
  df_list_genes[[drug]] <- data_list_genes
  
  
}

results_df <- bind_rows(results)

pad_vec <- function(x, n) {
  length(x) <- n
  x
}

flat <- list()

for (drug in names(df_list_genes)) {
  inner <- df_list_genes[[drug]]
  for (nm in names(inner)) {
    flat[[nm]] <- inner[[nm]]
  }
}

max_len <- max(lengths(flat))

flat_padded <- lapply(flat, pad_vec, n = max_len)

df <- as.data.frame(flat_padded, stringsAsFactors = FALSE)


write.xlsx(df, "./data/result_gene_hits_change/olivieri_hits_original_refined_approved_symbol.xlsx")

