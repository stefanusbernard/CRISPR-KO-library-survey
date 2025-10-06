# sgRNA log2fc file sources

1. Avana library (AvanalogfoldChange.csv, DepMap 25Q2) --> https://depmap.org/portal/data_page/?tab=allData

2. TKOv3 library log-fold-change data obtained from  --> https://academic.oup.com/g3journal/article/7/8/2719/6031511

3. TKOv3 library read count data obtained from Identifying chemogenetic interactions from CRISPR screens with drugZ (RPE1 cells drugZ; matrix-reads-by-gRNA-RPE1-drugZ.txt) --> https://link.springer.com/article/10.1186/s13073-019-0665-3#availability-of-data-and-materials & https://figshare.com/articles/dataset/Readcounts/8799215?file=16170896 
    - List of read counts:
        - RPE1_T0	
        - RPE1_T3A_CTRL_LOPRIMER	
        - RPE1_T3A_CTRL	
        - RPE1_T3A_ENTI_LOPRIMER	
        - RPE1_T3A_ENTI	
        - RPE1_T3A_GEMC	
        - RPE1_T3A_VINC	
        - RPE1_T3B_CTRL	
        - RPE1_T3B_ENTI	
        - RPE1_T3B_GEMC	
        - RPE1_T3B_VINC
    - Section A general-use algorithm for drug-gene interactions
        - We further conducted an independent screen of hTERT-immortalized RPE1 epithelial cells to determine genetic modifiers of the microtubule stabilizing agent vincristine, they found ABCC1 a known marker for clinical resistance to vincristine is the top synthetic hit in their screen
        - We screened hTERT-RPE1 cells with gemcitabine, a pyrimidine nucleoside analog and drugZ reveals a synthetic lethal interaction with DTYMK
        - Entinostat was in the Biorxiv paper of DrugZ (https://www.biorxiv.org/content/10.1101/232736v2.full.pdf) and removed in the published version (https://link.springer.com/article/10.1186/s13073-019-0665-3) 

4. TKOv3 library (Olivieri et al, 2020) --> https://data.mendeley.com/datasets/gfcn2wmrpf/1
    - Dataset_S2_readcounts.txt (we only take Dataset_S2 because this read count data for drugs screened with TKOv3)
    - Screen_Library.csv (not all treatment/replicates were performed in TKOv3!)
    - TKOv3 containing 71090 guides targeting 18056 genes (Olivieri et al, 2020)

5. Brunello library (Samson et al, 2018) --> https://www.nature.com/articles/s41467-018-07901-8#article-info
    - The Brunello library comprises 77,441 sgRNAs, an average of 4 sgRNAs per gene, and 1000 non-targeting control sgRNAs. 
    - We conducted genome-wide negative selection (dropout) screens in A375 melanoma cells that were first engineered to express Cas9. 
