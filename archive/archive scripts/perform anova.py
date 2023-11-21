# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 15:49:53 2023

@author: robbi
"""

import alphastats 
import pandas as pd
import matplotlib.pyplot as plt

path = "G:/My Drive/Data/data/spikein data/"
meta = f"{path}meta.xlsx"
group = 'Treatment'

def preform_anova(path, meta, group):
    loader_total = alphastats.GenericLoader(f"{path}total_imputed.csv", 
                                      intensity_column = '[Run]',
                                       index_column="Protein.Group")

    loader_light = alphastats.GenericLoader(f"{path}light_imputed.csv", 
                                      intensity_column = '[Run]',
                                       index_column="Protein.Group")

    loader_nsp = alphastats.GenericLoader(f"{path}nsp_imputed.csv", 
                                      intensity_column = '[Run]',
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

    total_anova = df_total.anova(column = group)

    nsp_anova = df_nsp.anova(column = group)

    light_anova  = df_light.anova(column = group)

    total_anova.to_csv(f'{path}total_anova.csv',sep=',')
    light_anova.to_csv(f'{path}light_anova.csv',sep=',')
    nsp_anova.to_csv(f'{path}nsp_anova.csv',sep=',')
preform_anova(path, meta, group)