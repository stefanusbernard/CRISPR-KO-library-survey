#!/usr/bin/bash

# original library TKOv3
python 05_drugz.py -i ../../data/sgrna_lfc_data/tkov3_data/matrix-reads-by-gRNA-RPE1-drugZ.txt -o GEMC_output_original.txt -r LacZ,luciferase,EGFR -c RPE1_T3A_CTRL,RPE1_T3B_CTRL -x RPE1_T3A_GEMC,RPE1_T3B_GEMC

# refined library TKOv3
python 05_drugz.py -i ../../data/sgrna_lfc_data/tkov3_data/refined_matrix-reads-by-gRNA-RPE1-drugZ.txt -o GEMC_output_refined.txt -r LacZ,luciferase,EGFR -c RPE1_T3A_CTRL,RPE1_T3B_CTRL -x RPE1_T3A_GEMC,RPE1_T3B_GEMC
