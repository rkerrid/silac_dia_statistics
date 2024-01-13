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

# # path = 'G:/My Drive/Data/data/poc4/stats tests/normalize first/protein intensities/'
# # normalize_samples.main(path, group,'href')
# path = 'G:/My Drive/Data/data/poc4/H/normalized/'
# adapted_imputation.process_intensities(path, plot_imputation=True)
# path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed/'
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
path = 'G:/My Drive/Data/data/poc4/H/normalized/imputed/ttest results/'
pois = ['FTH1', 'FTL', 'TFRC', 'ACO1', 'IREB2','SLC11A2','CYBRD1', 'ACO2', 'BRD4', 'CDK4', 'CDK6','EIF4E','EIF4G1','SALL4', 'ZNF827','ZNF692','ZFP91','DTWD1', 'RNF166', 'FAM83F','EIF4EBP1']
loop_and_plot_results(path, pois, interactive=False, uniprot=False)