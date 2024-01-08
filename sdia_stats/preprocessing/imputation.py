# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 04:54:36 2023

@author: robbi
"""



from sdia_stats.utils.manage_directories import create_directory

from utils import manage_directories

import pandas as pd
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt

# from icecream import ic

 # import metadata
 

from icecream import ic
ic.enable()
 # import metadata
 

def import_meta(path):
     metadata = pd.read_csv(f"{path}meta.csv")
     return metadata
 
def get_dataframes(path, quantification):
    total = pd.read_csv(f"{path}protein intensities/total_{quantification}.csv", sep=',')
    light = pd.read_csv(f"{path}protein intensities/light_{quantification}.csv", sep=',')
    nsp = pd.read_csv(f"{path}protein intensities/nsp_{quantification}.csv", sep=',')
    return total, light, nsp
     
def replace_values(df):
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.iloc[:,1:] = df.iloc[:,1:].apply(lambda x: np.where(x < 0.001, np.nan, x))
    return df
 

def filter_for_valid_values(df, metadata): # problem with this function or method, work around for loss of proteins that are in the dataframe is setting to true, need to fix this
     df['keep_row'] = False

     for group in metadata['Treatment'].unique():
         sub_meta = metadata[metadata['Treatment'] == group]
         cols = sub_meta['Sample'].tolist()
def filter_for_valid_values(df, metadata): 
     df['keep_row'] = False
   
     for group in metadata['Treatment'].unique():
         sub_meta = metadata[metadata['Treatment'] == group]
         cols = sub_meta['Sample'].tolist()
         ic(cols)
         # Check that we have at least 2 columns to compare, if not continue to next group
         if len(cols) < 2:
             continue
 
         # Calculate the sum of valid (not NaN) values for the columns in cols
         valid_sum = df[cols].notna().sum(axis=1)

         # Update 'keep_row' to True where at least 2 of the columns are not NaN
         df.loc[valid_sum >= 2, 'keep_row'] = True
     df = df[df['keep_row']].copy()

     df = df[df['keep_row']]
     df.drop('keep_row', axis=1, inplace=True)
     return df     
 
# imputation
def create_distribution(data):
    mu, std = stats.norm.fit(data)
    mu = mu - 1.8 * std
    std = 0.25*std
   
    return mu, std
    
def impute(x, mu, std, imputed_only):
    
    if imputed_only:
        if np.isnan(x):
            return np.random.normal(mu, std)
        else: 
            return np.nan
        
    else:
        if np.isnan(x):
            return np.random.normal(mu, std)
        else:
            return x

def perform_imputation(df):
    cols = df.columns.values.tolist()[1:]
    imputed_values = pd.DataFrame(columns=df.columns.values.tolist())
    for condition in cols:
        data = df[condition].dropna()
        
        mu, std = create_distribution(data)
        
        imputed_values[condition] = df[condition].apply(lambda x: impute(x, mu, std, True))
        df[condition] = df[condition].apply(lambda x: impute(x, mu, std, False))
    
    return  df, imputed_values

def plot_histogram(df, imputed_values, title):
  # Plot histograms for each column in df1 and df2 on the same plot with different colors
    # Create a single figure with 2 rows and 3 columns
  
    # Plot histograms for each column in df and imputed_values
    for  col in df.columns.values.tolist()[1:]:
        imputed_data = imputed_values[col].dropna()
        original_data = df[col].dropna()
        
        plt.hist(original_data, bins=20, alpha=0.5, label='original data', color='blue')
        plt.hist(imputed_data, bins=20, alpha=0.5, label='imputed values', color='green')
        plt.title(f'{title} {col} Histogram')
        plt.xlabel('Log2 intensity')
        plt.ylabel('Frequency')
        # plt.legend()

  

        # Display the figure
        plt.show()

def subset_data(df,  metadata):
    # need to map the cols to keep witht he metadata and drop/save cols that meet re requirements

    relevant_samples = metadata['Sample'].values.tolist()
 
    columns_to_keep = ['Protein.Group'] + relevant_samples

    df = df[columns_to_keep]
    
    return df

def subset_metadata(metadata, subset):
    filtered_metadata = metadata[metadata['Treatment'].isin(subset)]
    return filtered_metadata


def process_intensities(path, subset=[], quantification='href', plot_imputation=False):
    
    metadata = import_meta(path)
    # groups = metadata[metadata_sample_group].unique()
    
    total, light, nsp = get_dataframes(path,quantification)
    if len(subset) > 0:
        metadata = subset_metadata(metadata, subset)
    total = subset_data(total,metadata)
    light = subset_data(light, metadata)
    nsp = subset_data(nsp, metadata)

def process_intensities(path, subset, quantification='href', plot_imputation=False):
    metadata = import_meta(path)
    # groups = metadata[metadata_sample_group].unique()
    
    total, light, nsp = get_dataframes(path,'href')
    print('imported')
    print(total[total['Protein.Group']=='P06730-EIF4E'])
    metadata = subset_metadata(metadata, subset)
    total = subset_data(total,metadata)
    light = subset_data(light, metadata)
    nsp = subset_data(nsp, metadata)
    print('subseted')
    print(total[total['Protein.Group']=='P06730-EIF4E'])
    # replace NaN and inf values
    total = replace_values(total)
    light = replace_values(light)
    nsp = replace_values(nsp)
  
    print('replace')
    print(total[total['Protein.Group']=='P06730-EIF4E'])
    
    # filter for valid values
    total = filter_for_valid_values(total, metadata)
    light = filter_for_valid_values(light, metadata)
    nsp = filter_for_valid_values(nsp, metadata)
    
    print('valid values')
    print(total[total['Protein.Group']=='P06730-EIF4E'])
    # log transform
    total.iloc[:,1:] = np.log2(total.iloc[:,1:])
    light.iloc[:,1:] = np.log2(light.iloc[:,1:])
    nsp.iloc[:,1:] = np.log2(nsp.iloc[:,1:])
    
    
    # impute with gausian shift
    total_df, total_df_imputed = perform_imputation(total)
    nsp_df, nsp_df_imputed = perform_imputation(nsp)
    light_df, light_df_imputed = perform_imputation(light)
 
    
    if plot_imputation:
        plot_histogram(total_df, total_df_imputed, 'Total')
        plot_histogram(nsp_df, nsp_df_imputed, "NSP")
        plot_histogram(light_df, light_df_imputed, "Light")
        # base 2 exponentiation before saving
    total_df.iloc[:,1:] = 2**total_df.iloc[:,1:]
    nsp_df.iloc[:,1:] = 2**nsp_df.iloc[:,1:]
    light_df.iloc[:,1:] = 2**light_df.iloc[:,1:]
    
    create_directory(f"{path}", "imputed")
    metadata.to_csv(f"{path}/imputed/meta.csv", sep=',',index=False)
    nsp_df.to_csv(f"{path}imputed/nsp.csv", sep=',',index=False)
    light_df.to_csv(f"{path}imputed/light.csv", sep=',',index=False)
    total_df.to_csv(f"{path}imputed/total.csv", sep=',',index=False)



