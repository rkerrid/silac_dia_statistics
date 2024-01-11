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

def anova_clustermap(df, path, samples, title, pval=0.0005):
    filtered = df[df['ANOVA_pvalue']<pval]
    cols = ['Protein.Group']+samples
    filtered = filtered[cols]
    
    filtered.iloc[:,1:] = filtered.iloc[:,1:].apply(zscore, axis=1)
    filtered = filtered.set_index('Protein.Group')
    
    # Define the figure size, width and especially height to accommodate more y-axis labels
    height = len(filtered.index) * 0.3  # Adjust 0.3 to a smaller number if labels overlap, or larger for more space
    width = 10
    
    g = sns.clustermap(filtered, figsize=(width, height), method='average', metric='euclidean',  cmap='viridis')
    
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=-45, fontsize=6)
    g.savefig(f'{path}cluster_{title}_{pval}.pdf')
    plt.show()
   
def merge_data(df, anova):
    merged_df = pd.merge(df, anova, on='Protein.Group')
    return merged_df

def plot_clustermap(path):
    print(path)
    light_df = pd.read_csv(f'{path}light.csv',sep=',', index_col=False)
    nsp_df = pd.read_csv(f'{path}nsp.csv',sep=',', index_col=False)
    nsp_anova = pd.read_csv(f"{path}nsp_anova.csv", sep=',', index_col=False)
    light_anova = pd.read_csv(f"{path}light_anova.csv", sep=',', index_col=False)
    
    nsp_merged = merge_data(nsp_df, nsp_anova)
    light_merged = merge_data(light_df, light_anova)
    
    samples = light_df.columns.values.tolist()[1:]
    anova_clustermap(nsp_merged, path, samples, "nsp")
    anova_clustermap(light_merged, path, samples, "light")

if __name__ =='__main__':
    path = 'G:/My Drive/Data/data/poc4/H/imputed/normalized/'
    path = 'G:/My Drive/Data/data/poc4/H/imputed/normalized/for_anova/'
    path = 'G:/My Drive/Data/data/eif4g optimization/imputed/normalized/for_anova/'
    path = 'G:/My Drive/Data/data/example poc/H/imputed/normalized/anova_results/'
    
    # samples = ['control_I', 'control_II', 'control_III', 'DFO_I','DFO_II','DFO_III']
    # samples = ['08a_1','08a_2','08a_3','18a_1','18a_2','18a_3','28a_1','28a_2','28a_3']
    
    plot_clustermap(path)

