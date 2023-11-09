# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 21:37:49 2023

@author: robbi
"""

import pandas as pd
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt




def anova_clustermap(df, anova_df, title):
    merged_df = pd.merge(df, anova_df, on='Protein.Group')
    merged_df = merged_df.iloc[:,1:]
    filtered = merged_df[merged_df['ANOVA_pvalue']<0.00001]
    
    filtered = filtered[['Protein.Group', 'eIF4E+ 8h 1', 'eIF4E+ 8h 2', 'eIF4E- 8h 1', 'eIF4E- 8h 2', 'eIF4G1+ 8h 1', 'eIF4G1+ 8h 2', 'eIF4G1- 8h 1', 'eIF4G1- 8h 2', 'eIF4G2+ 8h 1', 'eIF4G2+ 8h 2', 'eIF4G2- 8h 1', 'eIF4G2- 8h 2', 'eIF4G3+ 8h 1', 'eIF4G3+ 8h 2', 'eIF4G3- 8h 1', 'eIF4G3- 8h 2']]
    filtered.iloc[:,1:] = filtered.iloc[:,1:].apply(zscore)
    
    filtered = filtered.set_index('Protein.Group')
    
    
    # Define the figure size, width and especially height to accommodate more y-axis labels
    height = len(filtered.index) * 0.3  # Adjust 0.3 to a smaller number if labels overlap, or larger for more space
    width = 10
    
    g = sns.clustermap(filtered, figsize=(width, height), method='average', metric='euclidean',  cmap='viridis')
    
    
    # sns.heatmap(total_filtered.iloc[:400,:], annot=False, cmap='viridis', linewidths=.5, vmin=-1,vmax=1)
    
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=-45, fontsize=6)
    g.savefig(f'cluster_{title}.pdf')
    # plt.show()
    
path = "G:/My Drive/Data/data/eIF4F optimization/imputed intensities 8hr/"

total_df = pd.read_csv(f'{path}total.csv',sep=',', index_col=False)
light_df = pd.read_csv(f'{path}light.csv',sep=',')
nsp_df = pd.read_csv(f'{path}nsp.csv',sep=',')

total_anova = pd.read_csv("C:/phd projects/silac_dia_statistics/eIF4F optimization/total_anova.csv", sep=',', index_col=False)
nsp_anova = pd.read_csv("C:/phd projects/silac_dia_statistics/eIF4F optimization/nsp_anova.csv", sep=',', index_col=False)
light_anova = pd.read_csv("C:/phd projects/silac_dia_statistics/eIF4F optimization/light_anova.csv", sep=',', index_col=False)

anova_clustermap(total_df, total_anova, "Total")
anova_clustermap(nsp_df, nsp_anova, "nsp")
anova_clustermap(light_df, light_anova, "light")