# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 08:22:44 2023

@author: rkerrid
"""


import pandas as pd
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic

def anova_clustermap(df, path, samples, title):
    ic(df)
    filtered = df[df['ANOVA_pvalue']<0.01]
    ic(filtered)
    cols = ['Protein.Group']+samples
    ic(cols)
    filtered = filtered[cols]
    
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
    # merged_df = merged_df.iloc[:,1:]
    ic(merged_df)
    return merged_df

def plot_clustermap(path, samples):
    light_df = pd.read_csv(f'{path}light.csv',sep=',', index_col=False)
    nsp_df = pd.read_csv(f'{path}nsp.csv',sep=',', index_col=False)
    ic(light_df)
    nsp_anova = pd.read_csv(f"{path}nsp_anova.csv", sep=',', index_col=False)
    light_anova = pd.read_csv(f"{path}light_anova.csv", sep=',', index_col=False)
    
    
    nsp_merged = merge_data(nsp_df, nsp_anova)
    light_merged = merge_data(light_df, light_anova)
    
    anova_clustermap(nsp_merged, path, samples, "nsp")
    anova_clustermap(light_merged, path, samples, "light")

if __name__ =='__main__':
    path = 'G:/My Drive/Data/data/poc4/H/imputed/normalized/'
    path = 'G:/My Drive/Data/data/poc4/H/imputed/normalized/for_anova/'
    path = 'G:/My Drive/Data/data/eif4g optimization/imputed/normalized/for_anova/'
    
    samples = ['control_I', 'control_II', 'control_III', 'DFO_I','DFO_II','DFO_III']
    samples = ['08a_1','08a_2','08a_3','18a_1','18a_2','18a_3','28a_1','28a_2','28a_3']
    
    plot_clustermap(path, samples)

