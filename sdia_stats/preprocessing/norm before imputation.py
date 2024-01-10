# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:15:50 2023

@author: robbi
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
from sdia_stats.utils.manage_directories import create_directory
from functools import reduce



# get control and treatment pairs for light data based on metadata
def import_dataframe(path, quantification, channel):
    df = pd.read_csv(f'{path}{channel}_{quantification}.csv', sep=',')
    return df

def import_pairs(path, control, sample, quantification, channel):
    
    light_df = import_dataframe(path, quantification, channel )

    # Selecting columns that start with the specified  for control
    selected_columns = [col for col in light_df.columns if control in col]
    control_light = light_df[['Protein.Group'] + selected_columns]
    
    # Selecting columns that start with the specified  for sample
    selected_columns = [col for col in light_df.columns if sample in col]
    sample_light = light_df[['Protein.Group'] + selected_columns]   
    
    return control_light, sample_light

def filter_rows(control, sample):
    # step 1: remove non valid values from both dfs
    control = control.replace([0, np.inf, -np.inf, np.nan], np.nan)
    sample = sample.replace([0, np.inf, -np.inf, np.nan], np.nan)
    
    # Step 2: Count valid values (non-NaN) in each row
    valid_count_control = control.iloc[:, 1:].notnull().sum(axis=1)
    valid_count_sample = sample.iloc[:, 1:].notnull().sum(axis=1)
    
    # Step 3: Filter rows with at least 2 valid values
    control_filtered = control[valid_count_control >= 2]
    sample_filtered = sample[valid_count_sample >= 2]
    
    # # Step 4: Keep only rows where protein is present in both DataFrames
    # Step 1: Identify common proteins
    common_proteins = pd.Series(list(set(control_filtered['Protein.Group']) & set(sample_filtered['Protein.Group'])))
    
    # Step 2: Filter both DataFrames
    control_common = control_filtered[control_filtered['Protein.Group'].isin(common_proteins)]
    sample_common = sample_filtered[sample_filtered['Protein.Group'].isin(common_proteins)]
    
    return control_common, sample_common


def make_distribution(control, sample, title):
    # Calculate median for each shared protein in each dataframe
    control_median = control.groupby('Protein.Group').median().reset_index()
    sample_median = sample.groupby('Protein.Group').median().reset_index()
    
    # Convert to log2 space 
    control_median_log2 = np.log2(control_median.iloc[:, 1:].values.flatten())
    control_median_log2 = control_median_log2[~np.isnan(control_median_log2)]  # Remove NaN values
    
    sample_median_log2 = np.log2(sample_median.iloc[:, 1:].values.flatten())
    sample_median_log2 = sample_median_log2[~np.isnan(sample_median_log2)]  # Remove NaN values
    
    # Plotting
 
    plt.hist(control_median_log2, bins=50, alpha=0.8, color='blue', label='Control')
    plt.hist(sample_median_log2, bins=50, alpha=0.2, color='green', label='Sample')
    
    # Calculate and plot median lines
    control_median_value = np.median(control_median_log2)
    sample_median_value = np.median(sample_median_log2)
    
    plt.axvline(control_median_value, color='blue', linestyle='dashed', linewidth=2, label='Control Median')
    plt.axvline(sample_median_value, color='green', linestyle='dashed', linewidth=2, label='Sample Median')
    
    plt.xlabel('Log2 intensity')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Log2 Transformed Values for light control and {title}')
    plt.legend()
    plt.show()
    
    # calculate median shift to normalize distributions
    normalization_factor = control_median_value - sample_median_value
    sample_median_log2 += normalization_factor
    
    # Plotting
 
    plt.hist(control_median_log2, bins=50, alpha=0.8, color='blue', label='Control')
    plt.hist(sample_median_log2, bins=50, alpha=0.2, color='green', label='Sample')
    
    # Calculate and plot median lines
    control_median_value = np.median(control_median_log2)
    sample_median_value = np.median(sample_median_log2)
    
    plt.axvline(control_median_value, color='blue', linestyle='dashed', linewidth=2, label='Control Median')
    plt.axvline(sample_median_value, color='green', linestyle='dashed', linewidth=2, label='Sample Median')
    
    plt.xlabel('Log2 intensity')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Log2 Transformed Values for light and {title} to extract normalization factor')
    plt.legend()
    plt.show()
    
    return normalization_factor
    
# def make_distribution_based_on_shared_proteins(sample, control):
#     control_final, sample_filtered = filter_rows(sample, control)

def log_and_drop_nan_rows(df):
    df_numeric_cols = df.columns[1:]
    df[df_numeric_cols] = df[df_numeric_cols].replace([0, np.inf, -np.inf], np.nan)
    
    # Now apply log2 transformation
    df[df_numeric_cols] = np.log2(df[df_numeric_cols])
    
    # Remove rows where any value in the numeric columns is NaN
    # df[df_numeric_cols] = df.dropna(subset=df_numeric_cols)
    
    return df
def plot_normalization(df_control, df_sample, normalization_factor, title, nsp=False, control=False):
    if nsp:
        channel = 'nsp'
    else:
        channel = 'light'
    # get value cols
    control_numeric_cols = df_control.columns[1:]
    sample_numeric_cols = df_sample.columns[1:]
    
    df_control[control_numeric_cols] = df_control[control_numeric_cols].replace([0, np.inf, -np.inf], np.nan)
    df_sample[sample_numeric_cols] = df_sample[sample_numeric_cols].replace([0, np.inf, -np.inf], np.nan)

    # Calculate median values
    control_median_val = df_control[control_numeric_cols].median().median()
    sample_median_val = df_sample[sample_numeric_cols].median().median()      
    
    # Plot histogram before normalization 
    plt.hist(df_control[control_numeric_cols].values.flatten(), bins=50, alpha=0.8, color='blue', label='control')
    plt.axvline(control_median_val, color='blue', linestyle='dashed', linewidth=2, label='control Median')
    plt.hist(df_sample[sample_numeric_cols].values.flatten(), bins=50, alpha=0.8, color='green', label='Sample')
    plt.axvline(sample_median_val, color='green', linestyle='dashed', linewidth=2, label='Sample Median')
    
    plt.xlabel('Log2 intensity')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Log2 Transformed Values for {channel} {title} before normalizing')
    plt.legend()
    plt.show()
    
    ic(normalization_factor)
    # apply normalization
    if not control:
        df_control[control_numeric_cols] += normalization_factor
    df_sample[sample_numeric_cols] += normalization_factor
    
    # realculate median values
    control_median_val = df_control[control_numeric_cols].median().median()
    sample_median_val = df_sample[sample_numeric_cols].median().median() 
    
    # Plot histogram after normalization light
    plt.hist(df_control[control_numeric_cols].values.flatten(), bins=50, alpha=0.8, color='blue', label='control')
    plt.axvline(control_median_val, color='blue', linestyle='dashed', linewidth=2, label='control Median')
    plt.hist(df_sample[sample_numeric_cols].values.flatten(), bins=50, alpha=0.8, color='green', label='Sample')
    plt.axvline(sample_median_val, color='green', linestyle='dashed', linewidth=2, label='Sample Median')
    
    plt.xlabel('Log2 intensity')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of Log2 Transformed Values for {channel} {title} after normalizing')
    plt.legend()
    plt.show()
    

    
    return df_control, df_sample

def normalize_proteomes(control_light, sample_light, control_nsp, sample_nsp, normalization_factor, title):
      
      # log2 and drop nan
      control_light_transformed = log_and_drop_nan_rows(control_light)
      sample_light_transformed = log_and_drop_nan_rows(sample_light)
      control_nsp_transformed = log_and_drop_nan_rows(control_nsp)
      sample_nsp_transformed = log_and_drop_nan_rows(sample_nsp)
      
      # get value cols
      control_light_numeric_cols = control_light_transformed.columns[1:]
      control_light_transformed[control_light_numeric_cols] = control_light_transformed[control_light_numeric_cols].replace([0, np.inf, -np.inf], np.nan)
      
      sample_light_transformed_numeric_cols = sample_light_transformed.columns[1:]
      sample_light_transformed[sample_light_transformed_numeric_cols] = sample_light_transformed[sample_light_transformed_numeric_cols].replace([0, np.inf, -np.inf], np.nan)
      
      control_nsp_transformed_numeric_cols = control_nsp_transformed.columns[1:]
      control_nsp_transformed[control_nsp_transformed_numeric_cols] = control_nsp_transformed[control_nsp_transformed_numeric_cols].replace([0, np.inf, -np.inf], np.nan)
      
      sample_nsp_transformed_numeric_cols = sample_nsp_transformed.columns[1:]
      sample_nsp_transformed[sample_nsp_transformed_numeric_cols] = sample_nsp_transformed[sample_nsp_transformed_numeric_cols].replace([0, np.inf, -np.inf], np.nan)
      
      # apply and plot nnormalization
      control_light, control_nsp = plot_normalization(control_light_transformed, sample_light_transformed, normalization_factor, title, nsp=False, control=True)
      sample_light, sample_nsp = plot_normalization(control_light_transformed, sample_light_transformed, normalization_factor, title, nsp=True)

      return sample_light, control_nsp, sample_nsp
   
def reverse_log(df):
    df.iloc[:,1:] = 2**df.iloc[:,1:]
    return df

def join_dfs(left, right):
    return pd.merge(left, right, on='Protein.Group', how='outer')


def main(path, group, quantification):
    metadata = pd.read_csv(f'{path}../meta.csv',sep=',')
    nsp_df_list = []
    light_df_list = []
    for key, value in group.items():
        control = key
        samples = value
        
        
        for treatment in samples:
            # import light control and sample dfs
            control_light, sample_light = import_pairs(path, control, treatment, quantification, 'light')
            # import nsp control and sample dfs
            control_nsp, sample_nsp = import_pairs(path, control, treatment, quantification, 'nsp')
            
            # filter for valid ros of atleast 2 proteins in each dataframe for a given protein
            control_light_filtered, sample_light_filtered = filter_rows(control_light, sample_light)
            
            # make distribution for light proteomes from control and sample and normalize    
            normalization_factor = make_distribution(control_light_filtered, sample_light_filtered, treatment)
            
            # use normalization factor to normalize both light and corrosponding nsp datasets
            sample_light_transformed, sample_nsp_transformed, control_nsp_transformed = normalize_proteomes(control_light, sample_light, control_nsp, sample_nsp, normalization_factor, treatment)
            
            # reverse log
            sample_light = reverse_log(sample_light_transformed)
            sample_nsp = reverse_log(sample_nsp_transformed)
            # control_nsp = reverse_log(control_nsp_transformed)
            # cotrol_light = reverse_log(control_light)
            
            # add dfs to respective lists
            light_df_list.append(sample_light)
            nsp_df_list.append(sample_nsp)
        
        # initialize df list with controll samples
        light_df = pd.read_csv(f'{path}light_{quantification}.csv', sep=',')
        # Selecting control cols
        selected_columns = [col for col in light_df.columns if 'control' in col]
        control_df = light_df[['Protein.Group'] + selected_columns]
        light_df_list.append(control_df) 
        
        #normalize nsp control
        nsp_df = pd.read_csv(f'{path}nsp_{quantification}.csv', sep=',')
        selected_columns = [col for col in nsp_df.columns if 'control' in col]
        nsp_df = nsp_df[['Protein.Group'] + selected_columns]
        nsp_df_list.append(control_df) 
        
        # Apply the function across all dataframes in the list
        nsp_complete_df = reduce(join_dfs, nsp_df_list)
        light_complete_df = reduce(join_dfs, light_df_list)
        
        new_path = path + '../'
        create_directory(f"{new_path}", "normalized")
        metadata.to_csv(f"{new_path}normalized/meta.csv", sep=',', index=False)
        light_complete_df.to_csv(f"{new_path}normalized/light.csv", sep=',',index=False)
        nsp_complete_df.to_csv(f"{new_path}normalized/nsp.csv", sep=',',index=False)

            

if __name__ == '__main__':
    path = 'G:/My Drive/Data/data/poc4/H/protein intensities/'
    group = {'control':['FAC','DFO']}
    main(path, group,  'href')

    
 

