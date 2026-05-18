import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# function to normalize sgRNA LFC by median of all sgRNA (ignoring the NA)
def normalize_column(col):
    median = np.nanmedian(col)
    centered = col - median
    mad = np.nanmedian(np.abs(centered))
    if mad == 0 or np.isnan(mad):
        return centered
    return centered / mad

# function to scaled LFC by the absolute average LFC value across cell lines for core essential genes by Traver Hart (2014) https://www.embopress.org/doi/full/10.15252/msb.20145216
def scale_essential(df, hart_data_dir):
    list_common_essential = pd.read_csv(hart_data_dir, header=None)
    list_common_essential = list(list_common_essential[0])
    list_common_essential

    # subset guides targeting essential genes
    ess_guides = df.loc[df.index.get_level_values('gene').isin(list_common_essential)]

    # compute average LFC across cell lines for each essential guides
    ess_avg = ess_guides.mean(axis = 1)

    # take absolute value and average across guides
    scale_factor = np.mean(np.abs(ess_avg))
    print('Scale factor: ', scale_factor)

    df_scaled = df / scale_factor
    
    return df_scaled, list_common_essential

def sanity_check_scale_essential(df, df_scaled, list_common_essential):
    # Sanity check for the scaling LFC method
    ess_before = df.loc[df.index.get_level_values('gene').isin(list_common_essential)].stack().values
    ess_after = df_scaled.loc[df_scaled.index.get_level_values('gene').isin(list_common_essential)].stack().values

    plt.boxplot([ess_before, ess_after], labels = ['Before', 'After'])
    plt.axhline(y=0, color='r', linestyle='--')
    plt.ylabel('logFC guides targeting common essential genes')
    plt.show()
    
# function to correct LFC for gene copy number alterations 

# function to average the row of the spacer
def mean_row(df):
    df['mean'] = df.mean(axis = 1)
    df = df['mean']
    return df
