# sgRNA log2fc file sources

1. Avana library (AvanalogfoldChange.csv, DepMap 25Q2) --> https://depmap.org/portal/data_page/?tab=allData

2. TKOv3 library (hap1.foldchange) --> https://academic.oup.com/g3journal/article/7/8/2719/6031511

3. TKOv3 library (RPE1 cells drugZ; matrix-reads-by-gRNA-RPE1-drugZ.txt) --> https://link.springer.com/article/10.1186/s13073-019-0665-3#availability-of-data-and-materials & https://figshare.com/articles/dataset/Readcounts/8799215?file=16170896 

4. TKOv3 library (Olivieri et al, 2020) --> https://data.mendeley.com/datasets/gfcn2wmrpf/1
    - Dataset_S2_readcounts.txt (we only take Dataset_S2 because this read count data for drugs screened with TKOv3)
    - Screen_Library.csv (not all treatment/replicates were performed in TKOv3!)
    - TKOv3 containing 71090 guides targeting 18056 genes (Olivieri et al, 2020)
