# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 10:13:59 2023

@author: robbi
"""

import alphastats 
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
from sdia_stats.utils.manage_directories import create_directory
from icecream import ic


def generate_alphastats_objects(path, meta):
    meta_file = pd.read_csv(meta, sep=',')
    intensity_cols = meta_file['Sample'].values.tolist()
    loader_total = alphastats.GenericLoader(f"{path}total.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

    loader_nsp = alphastats.GenericLoader(f"{path}nsp.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

    loader_light = alphastats.GenericLoader(f"{path}light.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

    df_total = alphastats.DataSet(
        loader = loader_total,
        metadata_path = meta,
        sample_column = 'Sample')

    df_light = alphastats.DataSet(
        loader = loader_light,
        metadata_path = meta,
        sample_column = 'Sample')

    df_nsp = alphastats.DataSet(
        loader = loader_nsp,
        metadata_path = meta,
        sample_column = 'Sample')

    # preprocess data, log2 transform
    df_nsp.preprocess(subset=True)
    df_total.preprocess(subset=True)
    df_light.preprocess(subset=True)
    
    return df_nsp, df_light, df_total

def save_results(path, ttest, file_name):
    create_directory(f'{path}', 'ttest results')
    ttest.to_csv(f"{path}/ttest results/{file_name}",sep=',',index=False)

def ttest(path, meta, groups):
    nsp, light, total = generate_alphastats_objects(path, meta)

    # Iterate over the dictionary of groups
    for comparison_name, (group1, group2) in groups.items():
        # ttest for total
        result_total = total.diff_expression_analysis(column='Treatment', group1=group1, group2=group2)
        save_results(path, result_total, f'{comparison_name}_total.csv')
        
        # ttest for light
        result_light = light.diff_expression_analysis(column='Treatment', group1=group1, group2=group2)
        save_results(path, result_light, f'{comparison_name}_light.csv')
        
        # ttest for nsp
        result_nsp = nsp.diff_expression_analysis(column='Treatment', group1=group1, group2=group2)
        save_results(path, result_nsp, f'{comparison_name}_nsp.csv')
    
if __name__ =='__main__':
    path = 'G:/My Drive/Data/data/poc1 statstest/H/imputed/'
    # path = 'G:/My Drive/Data/data/eIF4F pilot/imputed/'
    meta = f'{path}meta.csv'
    groups = {'FAC': ('FAC', 'control'),
              'both': ('both', 'control'),
              'gu': ('gu', 'control')}
    # intensity_cols = ['control_1', 'control_2', 'control_3', 'FAC_1', 'FAC_2', 'FAC_3','gu_1', 'gu_2', 'gu_3']
    ttest(path, meta, groups)
  

    