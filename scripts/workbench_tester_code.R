library(tidyverse)
library(VennDiagram)
library(ggvenn)
library(ggbreak)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.T2TCHM13v2.0)

source("./00_CRISPR-KO-library-survey-functions.R")
source('~/CRISPR-KO-GuideRefine/GuideRefine_functions.R')


comparisons_with_y <- label_positions %>%
  mutate(comparisons = list(list(c("multi-target guides", "perfect"))))

comparison_trend <- list(c("multi-target guides","perfect"))

x <- combined_data$mean[combined_data$alignment_bin == 1]  
stats <- boxplot.stats(x)$stats


combined_data %>%
  ggplot(aes(alignment, mean, fill = alignment)) +
  geom_boxplot(outlier.shape = NA, size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 1) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
  scale_x_discrete(labels = c("perfect" = "On-target \nsgRNA",
                              "multi-target guides" = "multi-target \nsgRNA")) +
  labs(title = "", 
       x = "", 
       y = "Averaged median sgRNA Log2FC") +
  theme(plot.margin = margin(10, 10, 25, 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 20, vjust = 0.5, margin = margin(r = 20)),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5), 
        line = element_line((size = 2), colour = 'black')) +
  facet_wrap(~library, ncol = 5) +
  stat_compare_means(comparisons = comparison_trend,
                     method = "wilcox.test",
                     method.args = list(alternative = "less"),
                     size = 5,
                     step.increase = 0.05,
                     tip.length = 0.01,
                     label = "p.format",
                     label.y = 1)








facet_labels = c("brunello" = "Brunello",
                 "toronto_v3" = "TKOv3",
                 "yusa" = "Yusa (Project Score)",
                 "avana" = "Avana",
                 "jacquere" = "Jacquere")

ggboxplot(combined_data, 
          x = "alignment", 
          y = "mean",
          fill = "alignment", 
          palette = "jco",
          facet.by = "library",
          outliers = FALSE) +
  stat_compare_means(comparisons = list(c("multi-target guides", "perfect")),
                     method = "wilcox.test",
                     method.args = list(alternative = "less"),
                     size = 5,
                     step.increase = 0.05,
                     tip.length = 0.01,
                     label = "p.format",
                     label.y = 1) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
  scale_fill_manual(
    values = c("perfect" = "#1f77b4", 
               "multi-target guides" = "#ff7f0e"),
    labels = c("On-target", "Multi-targeting")  # legend text
  ) +
  labs(title = "", 
       x = "", 
       y = "Averaged median sgRNA Log2FC",
       fill = "Alignment:") +
  theme(plot.margin = margin(10, 10, 25, 10),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 20, vjust = 0.5, margin = margin(r = 20)),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5), 
        line = element_line((size = 2), colour = 'black'),
        legend.text = element_text(size = 12)) +
  facet_wrap(~library, ncol = 5, labeller = labeller(library = facet_labels))
  
  