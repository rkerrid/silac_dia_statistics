# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 04:54:36 2023

@author: robbi
"""


from manage_directories import create_directory
import pandas as pd
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt
from icecream import ic


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
 
def filter_for_valid_values(df, metadata, sample_group):
     for group in metadata[sample_group].unique():
         sub_meta = metadata[metadata[sample_group] == group]
         cols = sub_meta['Sample'].tolist()

         # Check that we have at least 2 columns to compare, if not continue to next group
         if len(cols) < 2:
             continue
 
         # Calculate the sum of valid (not NaN) values for the columns in cols
         valid_sum = df[cols].notna().sum(axis=1)

         # Update 'keep_row' to True where at least 2 of the columns are not NaN
         df.loc[valid_sum >= 2, 'keep_row'] = True
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
        ic(np.all(np.isfinite(data)))
        ic(data)
        mu, std = create_distribution(data)
        
        imputed_values[condition] = df[condition].apply(lambda x: impute(x, mu, std, True))
        df[condition] = df[condition].apply(lambda x: impute(x, mu, std, False))
    
    return  df, imputed_values

def plot_histogram(df, imputed_values):
  # Plot histograms for each column in df1 and df2 on the same plot with different colors
    # Create a single figure with 2 rows and 3 columns
  
    # Plot histograms for each column in df and imputed_values
    for  col in df.columns.values.tolist()[1:]:
        

        plt.hist(df[col], bins=20, alpha=0.5, label='original data', color='blue')
        plt.hist(imputed_values[col], bins=20, alpha=0.5, label='imputed values', color='green')
        plt.title(col + ' Histogram')
        plt.xlabel('Log2 intensity')
        plt.ylabel('Frequency')
        # plt.legend()

  

        # Display the figure
        plt.show()
 
def process_intensities(path, metadata_sample_group, quantification='href', plot_imputation=False):
    metadata = import_meta(path)
    # groups = metadata[metadata_sample_group].unique()
    total, light, nsp = get_dataframes(path,quantification)
    # replace NaN and inf values
    total = replace_values(total)
    light = replace_values(light)
    nsp = replace_values(nsp)
    
    # filter for valid values
    total = filter_for_valid_values(total, metadata, metadata_sample_group)
    light = filter_for_valid_values(light, metadata, metadata_sample_group)
    nsp = filter_for_valid_values(nsp, metadata, metadata_sample_group)
    
    # log transform
    total.iloc[:,1:] = np.log2(total.iloc[:,1:])
    light.iloc[:,1:] = np.log2(light.iloc[:,1:])
    nsp.iloc[:,1:] = np.log2(nsp.iloc[:,1:])
    
    
    # impute with gausian shift
    total_df, total_df_imputed = perform_imputation(total)
    nsp_df, nsp_df_imputed = perform_imputation(nsp)
    light_df, light_df_imputed = perform_imputation(light)
    
    if plot_imputation:
        plot_histogram(total_df, total_df_imputed)
        plot_histogram(nsp_df, nsp_df_imputed)
        plot_histogram(light_df, light_df_imputed)
        # base 2 exponentiation before saving
    total_df.iloc[:,1:] = 2**total_df.iloc[:,1:]
    nsp_df.iloc[:,1:] = 2**nsp_df.iloc[:,1:]
    light_df.iloc[:,1:] = 2**light_df.iloc[:,1:]
    
    create_directory(f"{path}", "imputed")
    nsp_df.to_csv(f"{path}/imputed/nsp.csv", sep=',')
    light_df.to_csv(f"{path}imputed/light.csv", sep=',')
    total_df.to_csv(f"{path}imputed/total.csv", sep=',')
    
# import meta and dataframes
path = 'G:/My Drive/Data/data/eIF4F pilot/'
metadata_sample_group = 'Treatment'
 
 
process_intensities(path, metadata_sample_group)
 
