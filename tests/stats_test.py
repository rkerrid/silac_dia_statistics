# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 16:22:42 2024

@author: robbi
"""
from sdia_stats.preprocessing import adapted_imputation 
from sdia_stats.preprocessing import normalize_samples
from sdia_stats.preprocessing import filter_contams_and_non_valid_rows
from sdia_stats.statistics import ttest
from sdia_stats.visualization.interpret_ttest_results import loop_and_plot_results
import pandas as pd




################################# srp


# # # filter contams and valid values
meta = 'G:/My Drive/Data/data/20240410 SRP AID/meta.csv'
# # normalize samples
path = 'G:/My Drive/Data/data/20240410 SRP AID/protein_groups/'

filter_contams_and_non_valid_rows.filter_protein_intensities(path, meta)

### normalize 
group = {'54-': ['54+'], 
          '68-':['68+'],
          '72-':['72+']
              }
path = 'G:/My Drive/Data/data/20240410 SRP AID/protein_groups_statistics/'
normalize_samples.main(path, group, meta)

### imputation
control_samples = list(group.keys())
light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path,control_samples, meta, plot_imputation=True)

### ttest
path = 'G:/My Drive/Data/data/20240410 SRP AID/protein_groups_statistics/imputed/'
# meta = f'{path}meta.csv'
# # Set the treatments you would like to compare using the t-test
groups = {
    "54": ('54+','54-'),
    "68": ('68+','68-'),
    "72": ('72+','72-')
    }

ttest.ttest(path, meta, groups)
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/protein_groups_filtered/imputed/ttest results/'
path = 'G:/My Drive/Data/data/20240410 SRP AID//protein_groups_statistics/imputed/ttest results/'
pois = ["SRP54", "SRP68", "SRP72"]
# pois = ["EIF4E",  "EIF4G1", "EIF4G2", "EIF4G3"]
loop_and_plot_results(path, pois, interactive=False, uniprot=False)



# # # # filter contams and valid values
# meta = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/'
# # # normalize samples
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/protein_groups/'
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/3D only/protein_groups/'

# filter_contams_and_non_valid_rows.filter_protein_intensities(path, meta)

# group = {'3D_08-': ['3D_08+'], 
#           '3D_28-':['3D_28+'],
#           'G3_08-':['G3_08+'],
#           'G2_28-':['G3_28+'],
#           'G2_08-':['G2_08+'],
#           'G2_28-':['G2_28+']
#               }
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/protein_groups_filtered/'
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/3D only/protein_groups_filtered/'
# normalize_samples.main(path, group)
# control_samples = list(group.keys())
# light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path,control_samples, plot_imputation=True)
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/protein_groups_filtered/imputed/'
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/3D only/protein_groups_filtered/imputed/'
# meta = f'{path}meta.csv'
# # # Set the treatments you would like to compare using the t-test
# groups = {
#     "3D08 vs control": ('3D_08+','3D_08-'),
#     "3D28 vs control": ('3D_28+','3D_28-'),
#     "G308 vs control": ('G3_08+','G3_08-'),
#     "G328 vs control":('G3_28+','G3_28-'),
#     "G208 vs control": ('G2_08+','G2_08-'),
#     "G228 vs control":('G2_28+','G2_28-')
#     }

# ttest.ttest(path, meta, groups)
# # path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/protein_groups_filtered/imputed/ttest results/'
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/3D only/protein_groups_filtered/imputed/ttest results/'
# pois = ["EIF4E", "EIF4E2", "EIF4E3", "EIF4G1", "EIF4G2", "EIF4G3", "EIF2S1", "EIF2A", "EIF5", "EIF5B", "EIF3A", "EIF3C", "EIF3D", "EIF4A1", "EIF4A2", "EIF4EP1", "EIF4EP2", "EEF2", "EIF2AK3", "EIF2AK4"]
# # pois = ["EIF4E",  "EIF4G1", "EIF4G2", "EIF4G3"]
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)



###############################################
# meta = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/'
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/protein_groups/'

# filter_contams_and_non_valid_rows.filter_protein_intensities(path, meta)

# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/protein_groups_filtered/'
# group = {'G1_08-': ['G1_08+'], 
#           'G1_28-':['G1_28+'],
#           'E_08-':['E_08+'],
#           'E_28-':['E_28+']          
#               }

# normalize_samples.main(path, group)

# control_samples = list(group.keys())
# light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path,control_samples, plot_imputation=True)

# meta = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/meta.csv'
# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/protein_groups_filtered/imputed/'
# groups = {
#     "G108 vs control": ('G1_08+','G1_08-'),
#     "G128 vs control": ('G1_28+','G1_28-'),
#     "E08 vs control": ('E_08+','E_08-'),
#     "E28 vs control":('E_28+','E_28-')
    
#     }

# ttest.ttest(path, meta, groups)

# path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/protein_groups_filtered/imputed/ttest results/'
# pois = ["EIF4E", "EIF4E2", "EIF4E3", "EIF4G1", "EIF4G2", "EIF4G3", "EIF2S1", "EIF2A", "EIF5", "EIF5B", "EIF3A", "EIF3C", "EIF3D", "EIF4A1", "EIF4A2", "EIF4EP1", "EIF4EP2", "EEF2", "EIF2AK3", "EIF2AK4"]
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)


########################################################################
meta = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/'
path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/protein_groups/'
path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/H but using only L and M/protein_groups/'
path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/H/protein_groups/'
meta = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/H/meta.csv'

# filter data
filter_contams_and_non_valid_rows.filter_protein_intensities(path, meta)

# Sample normalization
path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/H/protein_groups_statistics/'
group = {'control': ['FAC', 'DFO']          
              }
normalize_samples.main(path, group, meta)

# Imputation
path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/H/protein_groups_statistics/normalized/'
control_samples = list(group.keys())
light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path, control_samples, meta, plot_imputation=True)


# ttest
path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/protein_groups_statistics/imputed/'
groups = {
    "FAC vs control": ('FAC','control'),
    "DFO vs control": ('DFO','control')
    
    }

ttest.ttest(path, meta, groups)

path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/protein_groups_filtered/imputed/ttest results/'
pois = ['FTL','FTH1', 'TFRC']
loop_and_plot_results(path, pois, interactive=False, uniprot=False)







# # filter contams and valid values
# metadata_path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/'
# filter_contams_and_non_valid_rows.filter_protein_intensities(path, metadata_path, 'dlfq')
# # normalize samples
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/'
# normalize_samples.main(path, group)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/'
# subset_list = ['FAC', 'DFO', 'FACwARV', '4EGI', 'FACwGu']
# control_samples = list(group.keys())
# light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path, control_samples, plot_imputation=True)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/imputed/'
# meta = f'{path}meta.csv'
# # Set the treatments you would like to compare using the t-test
# groups = {
#     "FAC vs control": ('FAC','control'),
#     "DFO vs control": ('DFO', 'control'),












# path = 'G:/My Drive/Data/data/eIF4F optimization/protein intensities/'

# group = {'e-8h': ['e+8h'],
#           'g1-8h': ['g1+8h'],
#           'g2-8h': ['g2+8h'],
#           'g3-8h': ['g3+8h']
         
#               }


# # filter contams and valid values
# metadata_path = 'G:/My Drive/Data/data/eIF4F optimization/'
# filter_contams_and_non_valid_rows.filter_protein_intensities(path, metadata_path, 'href')
# # normalize samples
# path = 'G:/My Drive/Data/data/eIF4F optimization/protein_groups/'
# normalize_samples.main(path, group)
# path = 'G:/My Drive/Data/data/eIF4F optimization/protein_groups/normalized/'
# group = {'e-8h': ['e+8h'],
#           'g1-8h': ['g1+8h'],
#           'g2-8h': ['g2+8h'],
#           'g3-8h': ['g3+8h']
         
#               }
# control_samples = list(group.keys())
# light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path,control_samples, plot_imputation=True)
# path = 'G:/My Drive/Data/data/eIF4F optimization/protein_groups/normalized/imputed/'
# meta = f'{path}meta.csv'
# # # Set the treatments you would like to compare using the t-test
# groups = {
#     "e+ vs control": ('e+8h','e-8h'),
#     "g1+ vs control": ('g1+8h','g1-8h'),
#     "g2+ vs control": ('g2+8h','g2-8h'),
#     "g3+ vs control": ('g3+8h','g3-8h')
#     }

# ttest.ttest(path, meta, groups)
# path = 'G:/My Drive/Data/data/eIF4F optimization/protein_groups/normalized/imputed/ttest results/'
# pois = ["EIF4E", "EIF4E2", "EIF4E3", "EIF4G1", "EIF4G2", "EIF4G3", "EIF2S1", "EIF2A", "EIF5", "EIF5B", "EIF3A", "EIF3C", "EIF3D", "EIF4A1", "EIF4A2", "EIF4EP1", "EIF4EP2", "EEF2", "EIF2AK3", "EIF2AK4"]
# pois = ["EIF4E",  "EIF4G1", "EIF4G2", "EIF4G3"]
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)

# path = 'G:/My Drive/Data/data/eif4g optimization/protein_groups/normalized/imputed/ttest results/'
# pois = ["EIF4E", "EIF4E2", "EIF4E3", "EIF4G1", "EIF4G2", "EIF4G3", "EIF2S1", "EIF2A", "EIF5", "EIF5B", "EIF3A", "EIF3C", "EIF3D", "EIF4A1", "EIF4A2", "EIF4EP1", "EIF4EP2", "EEF2", "EIF2AK3", "EIF2AK4"]
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)

# # pois = ['EIF4EBP1']
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H/protein_groups/normalized/imputed/ttest results/'
# pois = ['FTH1', 'FTL', 'TFRC', 'ACO1', 'IREB2','SLC40A1','FBXL5 ', 'STEAP3', 'HEPH']
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)

# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/new h filter test/protein intensities/'
# # filter contams and valid values
# metadata_path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/new h filter test/'
# filter_contams_and_non_valid_rows.filter_protein_intensities(path, metadata_path, 'href')
# # normalize samples
# group = {'control': ['FAC', 'DFO']}

# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/new h filter test/protein_groups/'
# normalize_samples.main(path, group)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/new h filter test/protein_groups/normalized/'
# control_samples = list(group.keys())
# subset_list = ['control','FAC', 'DFO']
# light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path, control_samples, subset_list, plot_imputation=True)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/new h filter test/protein_groups/normalized/imputed/'
# meta = f'{path}meta.csv'
# # Set the treatments you would like to compare using the t-test
# groups = {
#     "FAC vs control": ('FAC','control'),
#     "DFO vs control": ('DFO', 'control')}


# ttest.ttest(path, meta, groups)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/new h filter test/protein_groups/normalized/imputed/ttest results/'
# pois = ['FTH1', 'FTL', 'TFRC', 'ACO1', 'IREB2','SLC40A1','FBXL5 ', 'STEAP3', 'HEPH']
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)

# group = {'control': ['FAC', 'DFO', 'FACwARV', '4EGI', 'FACwGu']}

# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein intensities/'
# # filter contams and valid values
# metadata_path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/'
# filter_contams_and_non_valid_rows.filter_protein_intensities(path, metadata_path, 'dlfq')
# # normalize samples
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/'
# normalize_samples.main(path, group)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/'
# subset_list = ['FAC', 'DFO', 'FACwARV', '4EGI', 'FACwGu']
# control_samples = list(group.keys())
# light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = adapted_imputation.process_intensities(path, control_samples, plot_imputation=True)
# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/imputed/'
# meta = f'{path}meta.csv'
# # Set the treatments you would like to compare using the t-test
# groups = {
#     "FAC vs control": ('FAC','control'),
#     "DFO vs control": ('DFO', 'control'),
#     "4EGI-1 vs control":('4EGI','control'),
#     "FAC with ARV825 vs control": ('FACwARV','control')}

# ttest.ttest(path, meta, groups)


# path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/protein_groups/imputed/ttest results/'
# pois = ['FTH1', 'FTL', 'TFRC', 'ACO1', 'IREB2','SLC40A1','FBXL5 ', 'STEAP3', 'HEPH']
# loop_and_plot_results(path, pois, interactive=False, uniprot=False)



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


