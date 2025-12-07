library(tidyverse)
library(ggpubr)
library(BSgenome.Hsapiens.NCBI.T2TCHM13v2.0)
library(HGNChelper)

source('~/CRISPR-KO-library-survey/scripts/00_CRISPR-KO-library-survey-functions.R')


terms <- c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "LacZ", "luciferase", "NO_SITE", "ONE_INTERGENIC_SITE")
# lacZ and luciferase are the control in Toronto V3 library

terms_for_filtering <- paste(terms, collapse = "|")

brunello <- 'broadgpp-brunello-library-contents'
toronto <- 'tkov3_guide_sequence'
avana <- 'avana_library'
yusa <- 'yusa_hcrispr_ko_grnas'
jacquere <- 'Jacquere_PerGuideAnnotations_Quota4'

# import sgRNA refining report

avana_report <- import_and_process_report(paste0('~/CRISPR-KO-GuideRefine/output_cleaning/T2T-CHM13/', avana, '_full_report.xlsx'))
brunello_report <- import_and_process_report(paste0('~/CRISPR-KO-GuideRefine/output_cleaning/T2T-CHM13/', brunello, '_full_report.xlsx'))
toronto_report <- import_and_process_report(paste0('~/CRISPR-KO-GuideRefine/output_cleaning/T2T-CHM13/', toronto, '_full_report.xlsx'))
yusa_report <- import_and_process_report(paste0('~/CRISPR-KO-GuideRefine/output_cleaning/T2T-CHM13/', yusa, '_full_report.xlsx'))
jacquere_report <- import_and_process_report(paste0('~/CRISPR-KO-GuideRefine/output_cleaning/T2T-CHM13/', jacquere, '_full_report.xlsx'))



sum(toronto_report$`sgRNA number`)
sum(toronto_report$`Additional sgRNA (Corrected)`)
sum(toronto_report$`actual total sgRNA`)


sum(toronto_report$`multi-target sgRNA`)
sum(toronto_report$`single mismatch sgRNA`)
sum(toronto_report$`PAM distal double mismatch`)

round(calculate_percentage(toronto_report, "single mismatch sgRNA", "sgRNA number"), 2)












