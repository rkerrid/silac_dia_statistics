# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 10:13:59 2023

@author: robbi
"""

import alphastats 
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
from manage_directories import create_directory

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