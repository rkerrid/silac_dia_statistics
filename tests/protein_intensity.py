# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:49:49 2024

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
    path = f'{path}/protein_groups/'
    intensity_cols = meta_file['Sample'].values.tolist()
    # loader_total = alphastats.GenericLoader(f"{path}total.csv", 
    #                                   intensity_column = intensity_cols,
    #                                     index_column="Protein.Group")

    loader_nsp = alphastats.GenericLoader(f"{path}nsp.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

    loader_light = alphastats.GenericLoader(f"{path}light.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

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

    # preprocess data, log2 transform
    df_nsp.preprocess(subset=True)
    # df_total.preprocess(subset=True)
    df_light.preprocess(subset=True)
    
    return df_nsp, df_light  #, df_total


# # # filter contams and valid values
meta = 'G:/My Drive/Data/data/20240410 SRP AID/meta.csv'
# # normalize samples
path = 'G:/My Drive/Data/data/20240410 SRP AID/'


# nsp, light = generate_alphastats_objects(path,meta)

protein54 = 'P61011-SRP54'
protein68 = 'Q9UHB9-SRP68'
protein72 = 'O76094-SRP72'
# nsp.plot_intensity(protein_id=protein, group='54+')

light = pd.read_csv(f'{path}/protein_groups/light.csv', sep=',')
light54 = light[light['Protein.Group']==protein54]
light68 = light[light['Protein.Group']==protein68]
light72 = light[light['Protein.Group']==protein72]