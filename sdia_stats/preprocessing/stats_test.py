# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 16:22:42 2024

@author: robbi
"""
from sdia_stats.preprocessing import adapted_imputation 
from sdia_stats.preprocessing import imputation 
from sdia_stats.preprocessing import normalize_samples
from sdia_stats.preprocessing import filter_contams_and_non_valid_rows
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

# group = {'control': ['FAC', 'DFO', 'ARV4', 'FAC_ARV', 'ARV6','4EG1', 'FAC_Gu', 'Gu2', 'Gu8','Len1', 'Len10']
#               }

# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H/protein intensities/'
# # filter contams and valid values
# metadata_path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H/'
# filter_contams_and_non_valid_rows.filter_protein_intensities(path, metadata_path, 'href')
# # normalize samples
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H/protein_groups/'
# normalize_samples.main(path, group)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H/protein_groups/normalized/'
# light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path, plot_imputation=True)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H/protein_groups/normalized/imputed/'
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
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H/protein_groups/normalized/imputed/ttest results/'
# pois = ['FTH1', 'FTL', 'TFRC', 'ACO1', 'IREB2','SLC11A2','CYBRD1', 'ACO2', 'BRD4', 'CDK4', 'CDK6','EIF4E','EIF4G1','SALL4', 'ZNF827','ZNF692','ZFP91','DTWD1', 'RNF166', 'FAM83F','EIF4EBP1']
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)

group = {'control': ['FAC', 'DFO', 'ARV4', 'FACwARV', 'ARV6','4EGI', 'FACwGu', 'Gu2', 'Gu8','len1', 'len10']
              }

path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein intensities/'
# filter contams and valid values
metadata_path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/'
filter_contams_and_non_valid_rows.filter_protein_intensities(path, metadata_path, 'dlfq')
# normalize samples
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/'
# normalize_samples.main(path, group)
path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/'
subset_list = ['4EGI_1', '4EGI_2', '4EGI_3', 'ARV4_1', 'ARV4_2', 'ARV4_3', 'ARV6_1', 'ARV6_2', 'ARV6_3', 'DFO_1', 'DFO_2', 'DFO_3', 'FAC_1', 'FAC_2', 'FAC_3', 'FACwARV_1', 'FACwARV_2', 'FACwARV_3', 'FACwGu1', 'FACwGu_2', 'FACwGu_3',  'Gu8_1', 'Gu8_2', 'Gu8_3', 'control_1', 'control_2', 'control_3', 'len10_1', 'len10_2', 'len10_3', 'len1_1', 'len1_2', 'len1_3']
subset_list = ['4EGI', 'ARV4', 'ARV6', 'DFO','FAC','FACwARV','FACwGu','Gu8','control','len1','len10']
light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path, subset=subset_list,plot_imputation=True)
path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/imputed/'
meta = f'{path}meta.csv'
# Set the treatments you would like to compare using the t-test
groups = {
    "FAC vs control": ('FAC','control'),
    "DFO vs control": ('DFO', 'control'),
    "ARV825 vs control": ('ARV6','control'),
    "Gu3340 vs control":('Gu8', 'control'),
    "4EGI-1 vs control":('4EGI','control'),
    "FAC with ARV825 vs control": ('FACwARV','control')}

ttest.ttest(path, meta, groups)
path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/imputed/ttest results/'
pois = ['FTH1', 'FTL', 'TFRC', 'ACO1', 'IREB2','SLC11A2','CYBRD1', 'ACO2', 'BRD4', 'CDK4', 'CDK6','EIF4E','EIF4G1','SALL4', 'ZNF827','ZNF692','ZFP91','DTWD1', 'RNF166', 'FAM83F','EIF4EBP1']
loop_and_plot_results(path, pois, interactive=False, uniprot=False)



# from sdia_stats.statistics import anova
# path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed/'
# # subset data
# light_df = pd.read_csv(f'{path}light.csv', sep=',')
# nsp_df = pd.read_csv(f'{path}nsp.csv', sep=',')
# # import meta
# meta = pd.read_csv(f'{path}meta.csv', sep =',')
# #subset_cols = ['Protein.Group','26a_1','26a_2','26a_3','26c_1','26c_2','26c_3']
# subset_cols = ['Protein.Group','DFO_I','DFO_II','DFO_III','FAC_I','FAC_II','FAC_III']
# title = 'fac v dfo new imputation'
# light_df = light_df[subset_cols]
# nsp_df = nsp_df[subset_cols]


# from sdia_stats.utils.manage_directories import create_directory

# # creat new directory for subsetted data to store anova results
# # create subdirectory for all anova
# create_directory(f'{path}', 'anova_results')
# path = f'{path}/anova_results/'

# # create subdirectory for this analysis
# create_directory(f'{path}', title)

# # save files to this subdirectory (including copy of metadata
# path = f'{path}/{title}/'
# light_df.to_csv(f'{path}/light.csv', index=False)
# nsp_df.to_csv(f'{path}/nsp.csv', index=False)
# meta.to_csv(f'{path}/meta.csv', index=False)

# group = 'Treatment'

# anova.preform_anova(path, group, subset_cols, title)
# from sdia_stats.visualization import clustermap_anova
# path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed//anova_results/fac v dfo new imputation/'
# clustermap_anova.plot_clustermap(path)


