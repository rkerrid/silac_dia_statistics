# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 16:22:42 2024

@author: robbi
"""
from sdia_stats.preprocessing import adapted_imputation 
from sdia_stats.preprocessing import imputation 
from sdia_stats.preprocessing import normalize_samples
from sdia_stats.statistics import ttest
from sdia_stats.visualization.interpret_ttest_results import loop_and_plot_results
import pandas as pd

# group = {'control': ['FAC', 'DFO', 'ARV4', 'FAC_ARV', 'ARV6','4EG1', 'FAC_Gu', 'Gu2', 'Gu8','Len1', 'Len10']
#              }

# # run imputation first
# path = 'G:/My Drive/Data/data/poc4/stats tests/imput first/protein intensities/'
# imputation.process_intensities(path, plot_imputation=True)
# path = 'G:/My Drive/Data/data/poc4/stats tests/imput first/protein intensities/imputed/'
# normalize_samples.main(path, group, 'href')
# path = 'G:/My Drive/Data/data/poc4/stats tests/imput first/protein intensities/normalized/'

# meta = f'{path}meta.csv'
# # Set the treatments you would like to compare using the t-test
# groups = {
#     "FAC vs control": ('FAC','control'),
#     "DFO vs control": ('DFO', 'control'),
#     "ARV825 vs control": ('ARV6','control'),
#     "Gu3340 vs control":('Gu8', 'control'),
#     "4EGI-1 vs control":('4EG1','control'),
#     "FAC with ARV825 vs control": ('FAC_ARV','control')}

# ttest.ttest(path, meta, groups)
# path = 'G:/My Drive/Data/data/poc4/stats tests/imput first/protein intensities/normalized/ttest results/'
# pois = ['FTH1', 'FTL', 'TFRC', 'BRD4', 'CDK4', 'CDK6','EIF4E','EIF4G1','SALL4', 'ZNF827','ZNF692','ZFP91','DTWD1', 'RNF166', 'FAM83F']
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)


# run sample norm first

group = {'control': ['FAC', 'DFO', 'ARV4', 'FAC_ARV', 'ARV6','4EG1', 'FAC_Gu', 'Gu2', 'Gu8','Len1', 'Len10']
              }

# path = 'G:/My Drive/Data/data/poc4/stats tests/normalize first/protein intensities/'
# normalize_samples.main(path, group,'href')
path = 'G:/My Drive/Data/data/poc4/H/normalized/'
adapted_imputation.process_intensities(path, plot_imputation=True)
path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed/'
meta = f'{path}meta.csv'
# Set the treatments you would like to compare using the t-test
groups = {
    "FAC vs control": ('FAC','control'),
    "DFO vs control": ('DFO', 'control'),
    "ARV825 vs control": ('ARV6','control'),
    "Gu3340 vs control":('Gu8', 'control'),
    "4EGI-1 vs control":('4EG1','control'),
    "FAC with ARV825 vs control": ('FAC_ARV','control')}

ttest.ttest(path, meta, groups)
path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed/ttest results/'
pois = ['FTH1', 'FTL', 'TFRC', 'ACO1', 'IREB2','SLC11A2','CYBRD1', 'ACO2', 'BRD4', 'CDK4', 'CDK6','EIF4E','EIF4G1','SALL4', 'ZNF827','ZNF692','ZFP91','DTWD1', 'RNF166', 'FAM83F','EIF4EBP1']
loop_and_plot_results(path, pois, interactive=False, uniprot=False)

from sdia_stats.statistics import anova
path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed/'
# subset data
light_df = pd.read_csv(f'{path}light.csv', sep=',')
nsp_df = pd.read_csv(f'{path}nsp.csv', sep=',')
# import meta
meta = pd.read_csv(f'{path}meta.csv', sep =',')
#subset_cols = ['Protein.Group','26a_1','26a_2','26a_3','26c_1','26c_2','26c_3']
subset_cols = ['Protein.Group','DFO_I','DFO_II','DFO_III','FAC_I','FAC_II','FAC_III']
title = 'fac v dfo new imputation'
light_df = light_df[subset_cols]
nsp_df = nsp_df[subset_cols]


from sdia_stats.utils.manage_directories import create_directory

# creat new directory for subsetted data to store anova results
# create subdirectory for all anova
create_directory(f'{path}', 'anova_results')
path = f'{path}/anova_results/'

# create subdirectory for this analysis
create_directory(f'{path}', title)

# save files to this subdirectory (including copy of metadata
path = f'{path}/{title}/'
light_df.to_csv(f'{path}/light.csv', index=False)
nsp_df.to_csv(f'{path}/nsp.csv', index=False)
meta.to_csv(f'{path}/meta.csv', index=False)

group = 'Treatment'

anova.preform_anova(path, group, subset_cols, title)
from sdia_stats.visualization import clustermap_anova
path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed//anova_results/fac v dfo new imputation/'
clustermap_anova.plot_clustermap(path)


