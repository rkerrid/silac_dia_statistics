# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 11:58:31 2024

@author: robbi
"""


from sdia_stats.utils.manage_directories import create_directory

import pandas as pd
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt

# from icecream import ic

 # import metadata
 

from icecream import ic
ic.disable()
 # import metadata
 

def import_meta(path):
     metadata = pd.read_csv(f"{path}meta.csv")
     return metadata
 
def get_dataframes(path, quantification):
    # total = pd.read_csv(f"{path}protein intensities/total_{quantification}.csv", sep=',')
    light = pd.read_csv(f"{path}light.csv", sep=',')
    nsp = pd.read_csv(f"{path}nsp.csv", sep=',')
    return  light, nsp
     
def replace_values(df):
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.iloc[:,1:] = df.iloc[:,1:].apply(lambda x: np.where(x < 0.001, np.nan, x))
    return df

         
def filter_for_valid_values(df, metadata): 
     df['keep_row'] = False
   
     for group in metadata['Treatment'].unique():
         sub_meta = metadata[metadata['Treatment'] == group]
         cols = sub_meta['Sample'].tolist()
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
    
def impute(x, mu, std): # at this point x is a row[col]
    
    # if imputed_only:
    if np.isnan(x):
        x = np.random.normal(mu, std)
        return x, True
    else:
        return x, False
        
    # else:
    #     if np.isnan(x):
    #         return np.random.normal(mu, std)
    #     else:
    #         return x

# def perform_imputation(df):
#     cols = df.columns.values.tolist()[1:]
#     imputed_values = pd.DataFrame(columns=df.columns.values.tolist())
#     for condition in cols:
#         data = df[condition].dropna()
        
#         global_mu, global_std = create_distribution(data)
        
#         # here need to check the row and apply the imputation col wise but row at a time conditional on row wise min non NaN mean val
#         imputed_values[condition] = df[condition].apply(lambda x: impute(x, global_mu, global_std, True))
#         df[condition] = df[condition].apply(lambda x: impute(x, global_std, global_std, False))
    
#     return  df, imputed_values


def add_lowest_sample_mean_for_row(df, metadata):
    # Extract the sample columns (assuming they start with 'sample_')
    df['lowest_mean'] = np.nan
    ic(df)
    for group in metadata['Treatment'].unique():
        sub_meta = metadata[metadata['Treatment'] == group]
        cols = sub_meta['Sample'].tolist()
        # now itterate through sub df and if there are at least 2 valid vals, return the mean
        # if the current mean is not NaN or less than the new mean, the new lowest mean is the mean of that sample group for that protein
        for index, row in df.iterrows():
            new_mean = 0
            # Calculate the sum of valid (not NaN) values for the columns in cols
            valid_sum = row[cols].notna().sum()
            ic(valid_sum)
            if valid_sum >= 2:
                new_mean = row[cols].mean()
                ic(new_mean)
            else:
                new_mean = 10000000000000000000000000000000000000000000000000000000000000 # fix this
            if not np.isnan(row['lowest_mean']) or row['lowest_mean'] > new_mean:
                    ic(row['lowest_mean'])
                    row['lowest_mean'] = new_mean
                    ic(row['lowest_mean'])
    return df


# def perform_imputation(df, metadata):
#     cols = df.columns.values.tolist()[1:]
#     # annotate df with a column of the min non nan mean val for each row
#     df['was_imputed'] = False
#     ic(df)
#     df = add_lowest_sample_mean_for_row(df, metadata)
#     ic(df)
#     for condition in cols:
#         data = df[condition].dropna()
#         global_mu, global_std = create_distribution(data)
        
        
#         for index, row in df.iterrows():
#             # is imputed greater than mean for that row?
#             if row['lowest_mean'] < global_mu:
#                 adjusted_mu = row['lowest_mean'] - global_mu
#                 row[condition], was_imputed = impute(row[condition], adjusted_mu, global_std)
#             else:
#                 # Apply adjusted mu and std
#                 row[condition], was_imputed = impute(row[condition], global_mu, global_std)
#             df['was_imputed'] = was_imputed

#     imputed_df = df[df['was_imputed']==True]
#     ic(imputed_df)
#     ic(df)
#     return df, imputed_df

def perform_imputation(df, metadata):
    cols = df.columns.values.tolist()[1:]
    df['was_imputed'] = False
    df = add_lowest_sample_mean_for_row(df, metadata)
    
    for condition in cols:
        data = df[condition].dropna()
        global_mu, global_std = create_distribution(data)
        
        for index, row in df.iterrows():
            if row['lowest_mean'] < global_mu:
                adjusted_mu = row['lowest_mean'] - global_mu
                imputed_value, was_imputed = impute(row[condition], adjusted_mu, global_std)
            else:
                imputed_value, was_imputed = impute(row[condition], global_mu, global_std)
            
            df.at[index, condition] = imputed_value
            if was_imputed and df.at[index, 'was_imputed'] == True:
                df.at[index, 'was_imputed'] = True

    imputed_df = df[df['was_imputed']]
    return df, imputed_df






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
    if len(subset)>0:
        metadata = metadata[metadata['Treatment'].isin(subset)]
    return metadata




def process_intensities(path, subset = [], plot_imputation=False, quantification='href'):
    metadata = import_meta(path)
    # groups = metadata[metadata_sample_group].unique()
    
    light, nsp = get_dataframes(path, quantification)
    print('imported')
    metadata = import_meta(path)
    metadata = subset_metadata(metadata, subset)
    light = subset_data(light, metadata)
    nsp = subset_data(nsp, metadata)
    print('subseted')
    # replace NaN and inf values
    light = replace_values(light)
    nsp = replace_values(nsp)
  
    print('replace')
  
    
    # filter for valid values
    light = filter_for_valid_values(light, metadata)
    nsp = filter_for_valid_values(nsp, metadata)
    
    print('valid values')
  
    # log transform
    light.iloc[:,1:] = np.log2(light.iloc[:,1:])
    nsp.iloc[:,1:] = np.log2(nsp.iloc[:,1:])
    
    
    # impute with gausian shift
    nsp_df, nsp_df_imputed = perform_imputation(nsp, metadata)
    light_df, light_df_imputed = perform_imputation(light, metadata)
 
    
    if plot_imputation:
        plot_histogram(nsp_df, nsp_df_imputed, "NSP")
        plot_histogram(light_df, light_df_imputed, "Light")
        # base 2 exponentiation before saving
    nsp_df.iloc[:,1:] = 2**nsp_df.iloc[:,1:]
    light_df.iloc[:,1:] = 2**light_df.iloc[:,1:]
    
    create_directory(f"{path}", "imputed")
    metadata.to_csv(f"{path}/imputed/meta.csv", sep=',',index=False)
    nsp_df.to_csv(f"{path}imputed/nsp.csv", sep=',',index=False)
    light_df.to_csv(f"{path}imputed/light.csv", sep=',',index=False)


if __name__ == '__main__':
    path = 'G:/My Drive/Data/data/poc4/H/protein intensities/'
    path = 'G:/My Drive/Data/data/poc4/H/normalized/'
    process_intensities(path, plot_imputation=True)
