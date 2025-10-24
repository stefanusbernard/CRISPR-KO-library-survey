library(tidyverse)
library(VennDiagram)
library(ggvenn)
library(ggbreak)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.T2TCHM13v2.0)

source("./00_CRISPR-KO-library-survey-functions.R")
source('~/CRISPR-KO-GuideRefine/GuideRefine_functions.R')


combined_data


stat.test <- combined_data %>%
  group_by(library) %>%
  pairwise_wilcox_test(
    mean ~ alignment,     # paired comparison
    paired = TRUE,    # indicates repeated measures
    p.adjust.method = "bonferroni"
  )

combined_data <- combined_data %>%
  mutate(xpos = paste0(library, "_", alignment))


combined_data$xpos <- factor(combined_data$xpos, levels = unique(combined_data$xpos))

my_comp <- list(
  c("brunello_perfect","brunello_multi-target guides"),
  c("toronto_v3_perfect","toronto_v3_multi-target guides"),
  c("yusa_perfect","yusa_multi-target guides"),
  c("avana_perfect","avana_multi-target guides"),
  c("jacquere_perfect","jacquere_multi-target guides")
)

ggboxplot(
  combined_data, 
  x = "library", 
  y = "mean",
  color = "alignment",        # differentiate datasets by color
  palette = "jco",
  outlier.shape = NA) + 
  stat_compare_means(
    comparisons = my_comp,
    method = "wilcox.test",
    label = "p.signif",
    paired = FALSE,
    aes(group = alignment)
  )


