# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 11:26:02 2024

@author: robbi
"""

from sdia_stats.utils.manage_directories import create_directory
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from icecream import ic


def import_dataframe(path, quantification, channel):
    """
    Import a dataframe based on given parameters.
    """
    # df = pd.read_csv(f'{path}{channel}_{quantification}.csv', sep=',')
    df = pd.read_csv(f'{path}{channel}_{quantification}.csv', sep=',')
    return df

def save_filtered_data(light, nsp, metadata, path):
    new_path = f"{path}../"
    create_directory(new_path,'protein_groups')

    metadata.to_csv(f"{new_path}protein_groups/meta.csv", sep=',', index=False)
    light.to_csv(f"{new_path}protein_groups/light.csv", sep=',', index=False)
    nsp.to_csv(f"{new_path}protein_groups/nsp.csv", sep=',', index=False)
    
    
def log_transform(df):
    df.iloc[:, 1:] = np.log2(df.iloc[:, 1:])
    df.iloc[:, 1:] = df.iloc[:, 1:].replace([0, np.inf, -np.inf, np.nan], np.nan)
    return df

def reverse_log(df):
    df.iloc[:,1:] = 2**df.iloc[:,1:]
    return df
    
def filter_for_valid_values(df):
    """
    Filter rows in DataFrame based on valid values criteria.
    """

    # Creating a boolean mask for rows with at least 2 valid values
    valid_mask = df.iloc[:, 1:].notnull().sum(axis=1) >= 2

    # Filtering the DataFrame based on the mask
    filtered_df = df[valid_mask]

    # Getting the DataFrame with rows that were filtered out
    filtered_out_df = df[~valid_mask]

    return filtered_df, filtered_out_df

def filter_out_contams(df):
    
    # Create a contaminants mask based on the cont_ string and make sure all values are boolean
    contams_mask = df['Protein.Group'].str.contains('Cont_', case=False, na=False)
    if not all(isinstance(x, bool) for x in contams_mask):
        print("contams_mask contains non-boolean values:", contams_mask[~contams_mask.isin([True, False])])
    
    contams_df = df[contams_mask]  # Dataframe with only contaminants
    ic(contams_df)
    clean_df = df[~contams_mask]  # Dataframe without contaminants
    ic(clean_df)
    return clean_df, contams_df

def plot_histogram(df1, df2,title,label_df_1, label_df_2):
    # Flatten the DataFrames into single series
    data1 = df1.iloc[:, 1:].values.flatten()
    data2 = df2.iloc[:, 1:].values.flatten()

    # Remove NaN values
    data1 = data1[~np.isnan(data1)]
    data2 = data2[~np.isnan(data2)]

    # Plotting histograms
    plt.hist(data1, bins=50, alpha=0.8, color='blue', label=label_df_1)
    plt.hist(data2, bins=50, alpha=0.8, color='green', label=label_df_2)

    plt.xlabel('Log2 intensity')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Log2 Transformed Intensities After {title}')
    plt.legend()
    plt.show()
    
    
def filter_df(df):
    df, df_non_valid_rows = filter_for_valid_values(df)
    return df, df_non_valid_rows

def filter_protein_intensities(path, metadata_path, quantification):
    # import data
    metadata = pd.read_csv(f'{metadata_path}meta.csv', sep=',')
    df_light = import_dataframe(path, quantification, 'light')
    df_nsp = import_dataframe(path, quantification, 'nsp')
    # log transform data
    df_light = log_transform(df_light)
    df_nsp = log_transform(df_nsp)
    # filter for valid rows
    df_light, df_light_non_valid_rows = filter_df(df_light)
    plot_histogram(df_light, df_light_non_valid_rows,'filtering for valid vals light','remaining','non valid rows')
    df_nsp, df_nsp_non_valid_rows = filter_df(df_nsp)
    plot_histogram(df_nsp, df_nsp_non_valid_rows,'filtering for valid vals nsp','remaining','non valid rows')
    # filter for contams
    df_light, df_light_contams = filter_out_contams(df_light)
    plot_histogram(df_light, df_light_contams,'filtering out contams light','remaining','contams')
    nsp_df, df_nsp_contams = filter_out_contams(df_nsp)
    plot_histogram(nsp_df, df_nsp_contams,'filtering out contams nsp','remaining','contams')
    # reverse log
    light_df  = reverse_log(df_light)
    nsp_df  = reverse_log(nsp_df)
    
    save_filtered_data(light_df, nsp_df, metadata, path)
    
    
if __name__ == '__main__':
    path = 'G:/My Drive/Data/data/poc4/H/protein intensities/'
    metadata_path = 'G:/My Drive/Data/data/poc4/H/'
    metadata = pd.read_csv(f'{metadata_path}meta.csv', sep=',')
    
    filter_protein_intensities(path, metadata_path, 'href')
    
    
    
    
    
    
    