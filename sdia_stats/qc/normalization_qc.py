# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 13:42:10 2024

@author: robbi
"""

import alphastats 
import pandas as pd
import matplotlib as plt
import plotly.express as px


def generate_alphastats_objects(path, meta, normalization):
    meta_file = pd.read_csv(meta, sep=',')
    path = f'{path}/protein_groups/'
    intensity_cols = meta_file['Sample'].values.tolist()
    # loader_total = alphastats.GenericLoader(f"{path}total.csv", 
    #                                   intensity_column = intensity_cols,
    #                                     index_column="Protein.Group")

    loader_nsp = alphastats.GenericLoader(f"{path}nsp{normalization}.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

    loader_light = alphastats.GenericLoader(f"{path}light{normalization}.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")
    if normalization!= '_lfq':
        loader_h = alphastats.GenericLoader(f"{path}href{normalization}.csv", 
                                          intensity_column = intensity_cols,
                                            index_column="Protein.Group")
        df_h = alphastats.DataSet(
            loader = loader_h,
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
        
        return df_nsp, df_light, df_h
    else:
        # df_total = alphastats.DataSet(
        #     loader = loader_total,
        #     metadata_path = meta,
        #     sample_column = 'Sample')
    
        df_light = alphastats.DataSet(
            loader = loader_light,
            metadata_path = meta,
            sample_column = 'Sample')
    
        df_nsp = alphastats.DataSet(
            loader = loader_nsp,
            metadata_path = meta,
            sample_column = 'Sample')
    
        
    
        
        return df_nsp, df_light  #, df_total


def pca(path, meta, normalization):
    if normalization!= '_lfq':
        nsp, light, h = generate_alphastats_objects(path, meta, normalization)
        print('Light')
        light.plot_pca('Treatment').show()
        print('NSP')
        nsp.plot_pca('Treatment').show()
        print('Heavy channel')
        h.plot_pca('Treatment').show()
    else:
        nsp, light = generate_alphastats_objects(path, meta, normalization)
        print('Light')
        light.plot_pca('Treatment').show()
        print('NSP')
        nsp.plot_pca('Treatment').show()
        
        
def correlation(path, meta, normalization):
    custom_colorscale = [
    [0.0, "blue"],   # Start with blue
    [0.5, "white"],  # Transition to white
    [1.0, "red"]     # End with red
    ]
    if normalization!= '_lfq':
        nsp, light, h = generate_alphastats_objects(path, meta, normalization)
        print('Light')
        
        corr_matrix = light.mat.transpose().corr(method='spearman')
        plot = px.imshow(corr_matrix, color_continuous_scale='Agsunset')
        plot.show()
        print('NSP')
        corr_matrix = nsp.mat.transpose().corr(method='spearman')
        plot = px.imshow(corr_matrix, color_continuous_scale='Agsunset')
        plot.show()       
        print('Heavy channel')
        corr_matrix = h.mat.transpose().corr(method='spearman')
        plot = px.imshow(corr_matrix, color_continuous_scale='Agsunset')
        plot.show()
    else:
        nsp, light = generate_alphastats_objects(path, meta, normalization)
        print('Light')
        corr_matrix = light.mat.transpose().corr(method='spearman')
        plot = px.imshow(corr_matrix, color_continuous_scale='Agsunset')
        plot.show() 
        print('NSP')
        corr_matrix = nsp.mat.transpose().corr(method='spearman')
        plot = px.imshow(corr_matrix, color_continuous_scale='Agsunset')
        plot.show()        