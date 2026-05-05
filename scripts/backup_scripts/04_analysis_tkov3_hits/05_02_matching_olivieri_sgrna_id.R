source("./scripts/00_CRISPR-KO-library-survey-functions.R")
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(dplyr)

# Fetch sequence based on sgRNA coordinate from Olivieri paper

sgrna_coordinate <- read_tsv("./data/read_count_data/Dataset_S2_readcounts_olivieri_paper.txt") %>% select(sgRNA)

sgrna_parsed <- sgrna_coordinate %>%
  mutate(sgRNA_ID = sgRNA) %>%
  separate(sgRNA, into = c("chr", "rest"), sep = ":") %>%
  separate(rest, into = c("range", "gene", "strand"), sep = "_") %>%
  separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  filter(!chr %in% c("EGFP", "LacZ", "luciferase"))

# convert to granges

gr <- GRanges(
  seqnames = sgrna_parsed$chr,
  ranges   = IRanges(start = sgrna_parsed$start,
                     end = sgrna_parsed$end),
  strand   = sgrna_parsed$strand
)

# extract sequences from hg19

genome <- BSgenome.Hsapiens.UCSC.hg19
spacers <- getSeq(genome, gr)

# add sequence back to the table

sgrna_parsed$spacer <- as.character(spacers)
sgrna_parsed <- sgrna_parsed %>% select(sgRNA_ID, spacer)


# Import TKOv3 library 

terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase", "NO_SITE", "ONE_INTERGENIC_SITE")
# lacZ and luciferase are the control in Toronto V3 library
terms_for_filtering <- paste(terms, collapse = "|")

toronto_library <- import_sgrna_library("./data/library_data/original_library/tkov3_guide_sequence.tsv", terms_for_filtering)
refined_library <- import_sgrna_library("./data/library_data/refined_library/tkov3_guide_sequence_61K_refined.tsv", terms_for_filtering)

# Left join sgRNA_parsed to toronto_library by spacer and gene
toronto_library <- toronto_library %>%
  left_join(sgrna_parsed, join_by(spacer)) %>%
  drop_na()

# Left join sgRNA_parsed to refined library by spacer and gene (check this carefully as we corrected the symbol based on the alignment)
refined_library <- refined_library %>%
  left_join(sgrna_parsed, join_by(spacer)) %>%
  drop_na()


# Import Olivieri Screen readcount
olivieri_screen <- read_tsv("./data/read_count_data/Dataset_S2_readcounts_olivieri_paper.txt") %>% 
  dplyr::rename("sgRNA_ID" = "sgRNA") %>%
  select(-Gene)

# Import CRISPR-KO screen dataset
screen_info <- read_csv("./data/read_count_data/Screen_Library.csv") %>% 
  mutate(Screen = str_replace_all(Screen, "\\.", "-"),
         Screen = str_replace(Screen, "illudinS", "IlludinS"),
         Screen = str_replace(Screen, "HU-2", "HU-acute")) %>%
  dplyr::rename("drug" = "Screen") %>%
  filter(Library == "TKOv3")

olivieri_screen_data <- data.frame(screen_sample = colnames(olivieri_screen)) %>%
  mutate(olivieri_column_list = screen_sample,
         screen_sample = str_replace(screen_sample, "HU_acute", "HU-acute")) %>%
  separate(screen_sample, into = c("sample", "drug", "day", "replicate"), sep = "_") %>%
  filter(!sample == "sgRNA") %>%
  mutate(drug = str_replace(drug, "MLN", "MLN4924"),
         drug = str_replace(drug, "Illudin", "IlludinS")) %>%
  mutate(across(where(is.character), ~ replace_na(.x, "initial"))) %>%
  left_join(screen_info, join_by(drug)) %>%
  # remove the !drug == T0 if you would like to see the T0 count
  filter(!sample %in% c(screen_info |> pull(drug), "NT"),
         !drug == "T0") %>%
  mutate(Library = replace_na(Library, "TKOv3"),
         Nlib = replace_na(Nlib, 3))

olivieri_screen_data <- split(olivieri_screen_data, olivieri_screen_data$sample)


olivieri_screen_original <- olivieri_screen %>%
  left_join(toronto_library, join_by(sgRNA_ID)) %>%
  relocate(sgRNA, spacer, gene, .after = sgRNA_ID) %>%
  drop_na() %>%
  select(-c(sgRNA_ID, spacer))

write_tsv(olivieri_screen_original, "./data/read_count_data/Dataset_S2_readcounts_original.txt")

olivieri_screen_refined <- olivieri_screen %>%
  left_join(refined_library, join_by(sgRNA_ID)) %>%
  relocate(sgRNA, spacer, gene, .after = sgRNA_ID) %>%
  drop_na() %>%
  select(-c(sgRNA_ID, spacer))

write_tsv(olivieri_screen_refined, "./data/read_count_data/Dataset_S2_readcounts_refined.txt")




