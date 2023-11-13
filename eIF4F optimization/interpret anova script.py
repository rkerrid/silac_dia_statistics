# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 21:37:49 2023

@author: robbi
"""

import pandas as pd
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic

def anova_clustermap(df, path, title):
    ic(df)
    filtered = df[df['ANOVA_pvalue']<0.01]
    ic(filtered)
    filtered = filtered[['Protein.Group', 'Gu3340_1', 'Gu3340_2', 'Gu3340_3','control_1', 'control_2','control_3']]
    ic(filtered)
    filtered.iloc[:,1:] = filtered.iloc[:,1:].apply(zscore, axis=1)
    
    filtered = filtered.set_index('Protein.Group')
    ic(filtered)
    
    
    # Define the figure size, width and especially height to accommodate more y-axis labels
    height = len(filtered.index) * 0.3  # Adjust 0.3 to a smaller number if labels overlap, or larger for more space
    width = 10
    
    g = sns.clustermap(filtered, figsize=(width, height), method='average', metric='euclidean',  cmap='viridis')
    
    
    # sns.heatmap(total_filtered.iloc[:400,:], annot=False, cmap='viridis', linewidths=.5, vmin=-1,vmax=1)
    
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=-45, fontsize=6)
    g.savefig(f'{path}cluster_{title}.pdf')
    plt.show()
   
def merge_data(df, anova):
    merged_df = pd.merge(df, anova, on='Protein.Group')
    merged_df = merged_df.iloc[:,1:]
    
    return merged_df


path = "G:/My Drive/Data/data/spikein data/"

total_df = pd.read_csv(f'{path}total_imputed.csv',sep=',', index_col=False)
light_df = pd.read_csv(f'{path}light_imputed.csv',sep=',', index_col=False)
nsp_df = pd.read_csv(f'{path}nsp_imputed.csv',sep=',', index_col=False)

total_anova = pd.read_csv(f"{path}total_anova.csv", sep=',', index_col=False)
nsp_anova = pd.read_csv(f"{path}nsp_anova.csv", sep=',', index_col=False)
light_anova = pd.read_csv(f"{path}light_anova.csv", sep=',', index_col=False)

total_merged = merge_data(total_df, total_anova)
# nsp_merged = merge_data(nsp_df, nsp_anova)
light_merged = merge_data(light_df, light_anova)

anova_clustermap(total_merged, path, "Total")
anova_clustermap(nsp_merged, path, "nsp")
anova_clustermap(light_merged, path, "light")