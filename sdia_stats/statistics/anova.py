# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 08:23:54 2023

@author: rkerrid
"""

import alphastats 
import pandas as pd
import matplotlib.pyplot as plt
from sdia_stats.utils.manage_directories import create_directory


def preform_anova(path, group, subset_cols, title):
    meta = f'{path}meta.csv'
    meta_file = pd.read_csv(meta, sep=',')
    intensity_cols = subset_cols[1:]
    
    loader_light = alphastats.GenericLoader(f"{path}light.csv", 
                                      intensity_column = intensity_cols,
                                       index_column="Protein.Group")

    loader_nsp = alphastats.GenericLoader(f"{path}nsp.csv", 
                                      intensity_column = intensity_cols,
                                       index_column="Protein.Group")

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
    df_light.preprocess(subset=True)

    nsp_anova = df_nsp.anova(column = group)
    light_anova  = df_light.anova(column = group)

    light_anova.to_csv(f'{path}light_anova.csv',sep=',')
    nsp_anova.to_csv(f'{path}/nsp_anova.csv',sep=',')
