# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 12:41:56 2024

@author: robbi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
from sdia_stats.utils.manage_directories import create_directory
from icecream import ic


def import_dataframe(path, channel):
    """
    Import a dataframe based on given parameters.
    """
    # df = pd.read_csv(f'{path}{channel}_{quantification}.csv', sep=',')
    df = pd.read_csv(f'{path}{channel}.csv', sep=',')
    return df

def select_columns(df, keyword, path):
    """
    Select columns from a DataFrame based on a keyword.
    """
    metadata = pd.read_csv(f'{path}../meta.csv', sep=',')
    ic(metadata.columns.values.tolist())
    matching_samples = metadata[metadata['Treatment'] == keyword]['Sample'].tolist()

    # selected_columns = [col for col in df.columns if keyword in col]
    return df[['Protein.Group'] + matching_samples]

def filter_rows(df):
    """
    Filter rows in DataFrame based on valid values criteria.
    """
    df = df.replace([0, np.inf, -np.inf, np.nan], np.nan)
    valid_count = df.iloc[:, 1:].notnull().sum(axis=1)
    return df[valid_count >= 2]

def find_common_proteins(df1, df2):
    """
    Find common proteins in two DataFrames.
    """
    common_proteins = set(df1['Protein.Group']) & set(df2['Protein.Group'])
    df1 = df1[df1['Protein.Group'].isin(common_proteins)]
    df2 = df2[df2['Protein.Group'].isin(common_proteins)]
    return df1, df2

# Other functions like make_distribution, log_and_drop_nan_rows, plot_normalization, normalize_proteomes, reverse_log, join_dfs remain the same
def get_normalization_factor(df_control, df_sample, title):
    """
    Create and display histograms for log2 transformed values of control and sample data.
    """
    df_control_copy = df_control.copy(deep=True)
    df_sample_copy = df_sample.copy(deep=True)
    # Log2 transform and flatten the data, then remove NaN values
    def prepare_data_for_plot(df):
        data_log2 = np.log2(df.iloc[:, 1:].values.flatten())
        return data_log2[~np.isnan(data_log2)]

    control_data = prepare_data_for_plot(df_control_copy)
    sample_data = prepare_data_for_plot(df_sample_copy)

    # Plotting
    def plot_histogram(control_data, sample_data, state):
        plt.hist(control_data, bins=100, alpha=0.8, color='blue', label='control')
        plt.hist(sample_data, bins=100, alpha=0.8, color='green', label='sample')
        control_median_value = np.median(control_data)
        sample_median_value = np.median(sample_data)
        plt.axvline(control_median_value, color='blue', linestyle='dashed', linewidth=2, label=f'Median control')
        plt.axvline(sample_median_value, color='green', linestyle='dashed', linewidth=2, label=f'Median sample')
        
        plt.xlabel('Log2 intensity')
        plt.ylabel('Frequency')
        plt.title(f'Histogram of Log2 Transformed intensities for light control and light {title} {state}')
        plt.legend()
        plt.show()
        return control_median_value, sample_median_value

    control_median, sample_median = plot_histogram(control_data, sample_data, 'before normalization')

    normalization_factor = control_median - sample_median
    
    sample_data += normalization_factor
    
    # plot after applying normalization
    x,y = plot_histogram(control_data, sample_data, 'after normalization')
    
    return normalization_factor

def log_and_drop_nan_rows(df):
    """
    Apply log2 transformation to a DataFrame and drop rows with NaN values.
    """
    numeric_cols = df.columns[1:]
    df.loc[:,numeric_cols] = np.log2(df[numeric_cols].replace([0, np.inf, -np.inf], np.nan))
    return df.dropna(subset=numeric_cols)

def plot_normalization(df_control, df_sample, normalization_factor, title, channel='light'):
    """
    Plot histograms before and after applying normalization.
    """
    def plot(df_control, df_sample, title_suffix):
        df_control_numeric = df_control.iloc[:, 1:].values.flatten()
        df_sample_numeric = df_sample.iloc[:, 1:].values.flatten()
        
        plt.hist(df_control_numeric, bins=50, alpha=0.8, color='blue', label='Control')
        plt.axvline(np.median(df_control_numeric), color='blue', linestyle='dashed', linewidth=2, label='Median')
        plt.hist(df_sample_numeric, bins=50, alpha=0.8, color='green', label='Sample')
        plt.axvline(np.median(df_sample_numeric), color='green', linestyle='dashed', linewidth=2, label='Median')
        
        plt.xlabel('Log2 intensity')
        plt.ylabel('Frequency')
        plt.title(f'Histogram of {channel} {title} {title_suffix}')
        plt.legend()
        plt.show()

    plot(df_control, df_sample, 'before normalization')
    df_sample.iloc[:, 1:] += normalization_factor
    df_control.iloc[:, 1:] += normalization_factor
    plot(df_control, df_sample, 'after normalization')

def normalize_proteomes(df, normalization_factor, title, channel):
    """
    Normalize proteomes for a single dataset (either light or nsp) based on the calculated normalization factor.
    """
    def plot(df, df_unnorm, title):
        df_numeric = df.iloc[:, 1:].values.flatten()
        df_unnorm_numeric = df_unnorm.iloc[:, 1:].values.flatten()
        
        plt.hist(df_numeric, bins=50, alpha=0.8, color='green')
        plt.hist(df_unnorm_numeric, bins=50, alpha=0.8, color='blue')
        plt.axvline(np.median(df_numeric), color='green', linestyle='dashed', linewidth=2, label='Median norm')
  
        plt.axvline(np.median(df_unnorm_numeric), color='blue', linestyle='dashed', linewidth=2, label='Median unnorm')
        
        plt.xlabel('Log2 intensity')
        plt.ylabel('Frequency')
        plt.title(f'Histogram of {channel} {title} after applying normalization')
        plt.legend()
        plt.show()
        
    df = log_and_drop_nan_rows(df)
    df_unnorm = df.copy(deep=True)
    df.iloc[:, 1:] += normalization_factor
    
    plot(df, df_unnorm, title)
    
    return df


def reverse_log(df):
    """
    Reverse the log2 transformation.
    """
    df.iloc[:, 1:] = 2 ** df.iloc[:, 1:]
    return df

def join_dfs(left, right):
    """
    Join two DataFrames on the 'Protein.Group' column.
    """
    return pd.merge(left, right, on='Protein.Group', how='outer')

def derive_normalization_factors(path, group):
    normalization_factors = []
    control_df_list = []
    
    for control in group.keys():
        # Get control data that will not be normalized
        df_light = import_dataframe(path, 'light')
        df_light_control = select_columns(df_light, control, path)
        control_df_list.append(df_light_control)
    
        for treatment in group[control]:
            # Process for light data
            control_light, sample_light = import_pairs(path, control, treatment,  'light')
            control_light_common, sample_light_common = find_common_proteins(control_light, sample_light)
            
            # Calculate normalization factor for light data
            normalization_factor = get_normalization_factor(control_light_common, sample_light_common, f'{treatment}')
            result_tuple = (control, treatment, normalization_factor)
            normalization_factors.append(result_tuple)
    return normalization_factors, control_df_list

def apply_normalization(normalization_set, path, channel):
    control, treatment, normalization_factor = normalization_set
    df = import_dataframe(path,  channel)
    df_sample = select_columns(df, treatment, path)
    df_sample = normalize_proteomes(df_sample, normalization_factor, treatment, channel)
    df_sample = reverse_log(df_sample)
  
    return df_sample

def main(path, group, meta):
    metadata = pd.read_csv(f'{meta}', sep=',')
   
    normalization_factors, control_df_list = derive_normalization_factors(path, group)
    print(normalization_factors)
    nsp_df_list, light_df_list = [], []
    df_light = import_dataframe(path, 'light')
    df_nsp = import_dataframe(path, 'nsp')
    for control in group.keys():
        control_df = select_columns(df_light, control, path)
        control_nsp_df = select_columns(df_nsp, control, path)
        # need to normalize this control nsp set
        # for control in control_nsp_df.columns.values.tolist()[1:]:
        #    print(control)
        #    for control_set, in normalization_factors:
        light_df_list.append(control_df)
        nsp_df_list.append(control_nsp_df)
        
    for normalization_set in normalization_factors:
        light_norm = apply_normalization(normalization_set, path, 'light')
        nsp_norm = apply_normalization(normalization_set, path,  'nsp')
        

        nsp_df_list.append(nsp_norm)
        light_df_list.append(light_norm)

    # Combine and save the processed data
    save_combined_data(nsp_df_list, light_df_list, path, metadata)
    
# Reintroduced import_pairs function
def import_pairs(path, control, sample, channel):
    df = import_dataframe(path, channel)

    control_df = select_columns(df, control, path)
    sample_df = select_columns(df, sample, path)   

    return control_df, sample_df

# def select_control_data(path, quantification, channel):
#     df = pd.read_csv(f'{path}{channel}_{quantification}.csv', sep=',')
#     return select_columns(df, 'control', path)

def save_combined_data(nsp_df_list, light_df_list, path, metadata):
    nsp_complete_df = reduce(join_dfs, nsp_df_list)
    light_complete_df = reduce(join_dfs, light_df_list)
    
    # new_path = f"{path}../"
    create_directory(path,'normalized')
    

    # metadata.to_csv(f"{path}normalized/meta.csv", sep=',', index=False)
    light_complete_df.to_csv(f"{path}normalized/light.csv", sep=',', index=False)
    nsp_complete_df.to_csv(f"{path}normalized/nsp.csv", sep=',', index=False)

if __name__ == '__main__':
    path = 'G:/My Drive/Data/data/poc4/H/protein_groups/'
    group = {'control': ['FAC', 'DFO', 'ARV4', 'FAC_ARV', 'ARV6','4EG1', 'FAC_Gu', 'Gu2', 'Gu8']
             }
    main(path, group, 'href')