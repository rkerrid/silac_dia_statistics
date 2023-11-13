# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 15:53:17 2023

@author: robbi
"""
import pandas as pd
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic
import plotly.express as px
import numpy as np


def merge_data(df, ttest):
    merged_df = pd.merge(df, ttest, on='Protein.Group')
    # merged_df = merged_df.iloc[:,1:]
    
    return merged_df


def ttest_clustermap(df, path, title):
    ic(df)
    filtered = df[df['pval']<0.01]
    filtered = filtered.iloc[:, :-2]
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

def volcano(df):
    df['-log10(p_value)'] = -np.log10(df['pval'])  # Calculate the negative log of p-value
    ic(df)
    fig = px.scatter(df, x='log2fc', y='-log10(p_value)',
                     title='Volcano Plot',
                     labels={'log2fc': 'log2fc Size', '-log10(p_value)': '-log10(p-value)'})
    ic(fig)
    # Enhance the plot with additional customizations, such as highlighting significant points
    significant_threshold = -np.log10(0.05)  # For example, p < 0.05
    fig.add_hline(y=significant_threshold, line_dash="dash", line_color="red")
    
    fig.show()



# import imputed dfs
path = "G:/My Drive/Data/data/eIF4F optimization/imputed intensities 8hr/"
light = pd.read_csv(f'{path}light.csv', sep=',')
total = pd.read_csv(f'{path}total.csv', sep=',')
nsp = pd.read_csv(f'{path}nsp.csv', sep=',')
#import ttests
e_total = pd.read_csv(f'{path}e_total.csv', sep=',')
g1_total = pd.read_csv(f'{path}g1_total.csv', sep=',')
g2_total = pd.read_csv(f'{path}g2_total.csv', sep=',')
g3_total = pd.read_csv(f'{path}g3_total.csv', sep=',')

e_light = pd.read_csv(f'{path}e_light.csv', sep=',')
g1_light = pd.read_csv(f'{path}g1_light.csv', sep=',')
g2_light = pd.read_csv(f'{path}g2_light.csv', sep=',')
g3_light = pd.read_csv(f'{path}g3_light.csv', sep=',')

e_nsp = pd.read_csv(f'{path}e_nsp.csv', sep=',')
g1_nsp = pd.read_csv(f'{path}g1_nsp.csv', sep=',')
g2_nsp = pd.read_csv(f'{path}g2_nsp.csv', sep=',')
g3_nsp = pd.read_csv(f'{path}g3_nsp.csv', sep=',')
# merge ttests onto sub dfds based on cols

e_total = merge_data(total[['Protein.Group','eIF4E- 8h 1','eIF4E- 8h 2', 'eIF4E+ 8h 1', 'eIF4E+ 8h 2']], e_total)
g1_total = merge_data(total[['Protein.Group','eIF4G1- 8h 1','eIF4G1- 8h 2', 'eIF4G1+ 8h 1', 'eIF4G1+ 8h 2']], g1_total)
g2_total = merge_data(total[['Protein.Group','eIF4G2- 8h 1','eIF4G2- 8h 2', 'eIF4G2+ 8h 1', 'eIF4G2+ 8h 2']], g2_total)
g3_total = merge_data(total[['Protein.Group','eIF4G3- 8h 1','eIF4G3- 8h 2', 'eIF4G3+ 8h 1', 'eIF4G3+ 8h 2']], g3_total)


e_light = merge_data(light[['Protein.Group','eIF4E- 8h 1','eIF4E- 8h 2', 'eIF4E+ 8h 1', 'eIF4E+ 8h 2']], e_light)
g1_light = merge_data(light[['Protein.Group','eIF4G1- 8h 1','eIF4G1- 8h 2', 'eIF4G1+ 8h 1', 'eIF4G1+ 8h 2']], g1_light)
g2_light = merge_data(light[['Protein.Group','eIF4G2- 8h 1','eIF4G2- 8h 2', 'eIF4G2+ 8h 1', 'eIF4G2+ 8h 2']], g2_light)
g3_light = merge_data(light[['Protein.Group','eIF4G3- 8h 1','eIF4G3- 8h 2', 'eIF4G3+ 8h 1', 'eIF4G3+ 8h 2']], g3_light)


e_nsp = merge_data(nsp[['Protein.Group','eIF4E- 8h 1','eIF4E- 8h 2', 'eIF4E+ 8h 1', 'eIF4E+ 8h 2']], e_nsp)
g1_nsp = merge_data(nsp[['Protein.Group','eIF4G1- 8h 1','eIF4G1- 8h 2', 'eIF4G1+ 8h 1', 'eIF4G1+ 8h 2']], g1_nsp)
g2_nsp = merge_data(nsp[['Protein.Group','eIF4G2- 8h 1','eIF4G2- 8h 2', 'eIF4G2+ 8h 1', 'eIF4G2+ 8h 2']], g2_nsp)
g3_nsp = merge_data(nsp[['Protein.Group','eIF4G3- 8h 1','eIF4G3- 8h 2', 'eIF4G3+ 8h 1', 'eIF4G3+ 8h 2']], g3_nsp)

# clustermap
# ttest_clustermap(e_nsp, path, 'e_nsp')
# ttest_clustermap(g1_nsp, path, 'g1_nsp')
# ttest_clustermap(g2_nsp, path, 'g2_nsp')
# ttest_clustermap(g3_nsp, path, 'g3_nsp')

# ttest_clustermap(e_light, path, 'e_light')
# ttest_clustermap(g1_light, path, 'g1_light')
# ttest_clustermap(g2_light, path, 'g2_light')
# ttest_clustermap(g3_light, path, 'g3_light')

# ttest_clustermap(e_total, path, 'e_total')
# ttest_clustermap(g1_total, path, 'g1_total')
# ttest_clustermap(g2_total, path, 'g2_total')
# ttest_clustermap(g3_total, path, 'g3_total')


# volcano plot
volcano(e_nsp)


