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


def generate_alphastats_objects(path):
    loader_total = alphastats.GenericLoader(f"{path}total.csv", 
                                      intensity_column = '[Sample]',
                                       index_column="Protein.Group")

    loader_nsp = alphastats.GenericLoader(f"{path}nsp.csv", 
                                      intensity_column = '[Sample]',
                                       index_column="Protein.Group")

    loader_light = alphastats.GenericLoader(f"{path}light.csv", 
                                      intensity_column = '[Sample]',
                                       index_column="Protein.Group")

    df_total = alphastats.DataSet(
        loader = loader_total,
        metadata_path = f"{path}meta.csv",
        sample_column = 'Sample')

    df_light = alphastats.DataSet(
        loader = loader_light,
        metadata_path = f"{path}meta.csv",
        sample_column = 'Sample')

    df_nsp = alphastats.DataSet(
        loader = loader_nsp,
        metadata_path = f"{path}meta.csv",
        sample_column = 'Sample')

    # preprocess data, log2 transform
    df_nsp.preprocess(subset=True)
    df_total.preprocess(subset=True)
    df_light.preprocess(subset=True)
    
    return df_nsp, df_light, df_total

def save_results(path, ttest, file_name):
    create_directory(f'{path}', 'ttest results')
    ttest.to_csv(f"{path}/ttest results/{file_name}",sep=',',index=False)
    
def ttest(path):
    nsp, light, total = generate_alphastats_objects(path)
    
    
    # ttest for total
    e_total = total.diff_expression_analysis(column='Treatment', group1='e+8h', group2='e-8h')
    g1_total = total.diff_expression_analysis(column='Treatment', group1='g1+8h', group2='g1-8h')
    g2_total = total.diff_expression_analysis(column='Treatment', group1='g2+8h', group2='g2-8h')
    g3_total = total.diff_expression_analysis(column='Treatment', group1='g3+8h', group2='g3-8h')
    
    # ttest for light
    e_light = light.diff_expression_analysis(column='Treatment', group1='e+8h', group2='e-8h')
    g1_light = light.diff_expression_analysis(column='Treatment', group1='g1+8h', group2='g1-8h')
    g2_light = light.diff_expression_analysis(column='Treatment', group1='g2+8h', group2='g2-8h')
    g3_light = light.diff_expression_analysis(column='Treatment', group1='g3+8h', group2='g3-8h')
    
    # ttest for nsp
    e_nsp = nsp.diff_expression_analysis(column='Treatment', group1='e+8h', group2='e-8h')
    g1_nsp = nsp.diff_expression_analysis(column='Treatment', group1='g1+8h', group2='g1-8h')
    g2_nsp = nsp.diff_expression_analysis(column='Treatment', group1='g2+8h', group2='g2-8h')
    g3_nsp = nsp.diff_expression_analysis(column='Treatment', group1='g3+8h', group2='g3-8h')
    

    save_results(path, e_nsp, 'e_nsp.csv')
    save_results(path, g1_nsp, 'g1_nsp.csv')
    save_results(path, g2_nsp, 'g2_nsp.csv')
    save_results(path, g3_nsp, 'g3_nsp.csv')
    
    save_results(path, e_light, 'e_light.csv')
    save_results(path, g1_light, 'g1_light.csv')
    save_results(path, g2_light, 'g2_light.csv')
    save_results(path, g3_light, 'g3_light.csv')
    
    save_results(path, e_total, 'e_total.csv')
    save_results(path, g1_total, 'g1_total.csv')
    save_results(path, g2_total, 'g2_total.csv')
    save_results(path, g3_total, 'g3_total.csv')
    
    