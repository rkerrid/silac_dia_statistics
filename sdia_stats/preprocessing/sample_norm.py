# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 11:57:34 2023

@author: robbi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from icecream import ic


def load_and_prepare_data(file_path):
    """Load data from a CSV file and prepare it for analysis."""
    try:
        data = pd.read_csv(file_path, index_col='Protein.Group')
        data.replace([np.inf, -np.inf], np.nan, inplace=True)
        data.dropna(how='all', subset=data.columns, inplace=True)
        return data
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def plot_log2_distributions(data, title):
    """Plot log2 distributions of the data."""
    log2_data = np.log2(data)
    plt.figure(figsize=(20, 10))
    sns.boxplot(data=log2_data)
    plt.xticks(rotation=90)
    plt.ylabel('Log2 Expression')
    plt.title(title)
    plt.tight_layout()
    plt.show()
    return log2_data

def normalize(log2_data, medians):
    """Normalize the log2 data and convert back to original scale."""
    normalized_data = log2_data.subtract(medians)
    
    return normalized_data

# def plot_normalized_distributions(data, title):
#     """Plot distributions of the normalized data."""
#     data=np.log2(data)
#     plt.figure(figsize=(20, 10))
#     sns.boxplot(data=data)
#     plt.xticks(rotation=90)
#     plt.ylabel('Normalized Expression')
#     plt.title(title)
#     plt.tight_layout()
#     plt.show()

def plot_normalized_distributions(data, title):
    """Plot distributions of the normalized data."""
    log2_data = np.log2(data)
    plt.figure(figsize=(20, 10))
    sns.boxplot(data=data)
    plt.xticks(rotation=90)
    plt.ylabel('Normalized Expression')
    plt.title(title)
    plt.tight_layout()
    plt.show()
    

def save_to_csv(df, file_path):
    """Save DataFrame to a CSV file."""
    try:
        df.to_csv(file_path)
        print(f"Saved to {file_path}")
    except Exception as e:
        print(f"Error saving file: {e}")


def normalize_samples(path):
    # Directories and file paths
    base_dir = f'{path}imputed'
    file_path_light = os.path.join(base_dir, 'light.csv')
    file_path_nsp = os.path.join(base_dir, 'nsp.csv')
    save_dir = os.path.join(base_dir, 'normalized')
    
    # Load and process the first dataset
    data_light = load_and_prepare_data(file_path_light)
    if data_light is not None:
        # ic(data_light)
        log2_data_light = plot_log2_distributions(data_light, 'Log2 Distributions of Protein Expression Across Light Samples')
        medians_light = log2_data_light.median()
        normalized_light = normalize(log2_data_light, medians_light)
        # ic(normalized_light)
        normalized_light = np.power(2, normalized_light) # Inverse log2 transformation
        plot_normalized_distributions(normalized_light, 'Normalized Distributions of Protein Expression Across Light Samples')
    
    # Load and process the second dataset
    data_nsp = load_and_prepare_data(file_path_nsp)
    if data_nsp is not None:
        ic(data_nsp)
        log2_data_nsp = plot_log2_distributions(data_nsp, 'Log2 Distributions of Protein Expression Across NSP Samples')
        ic(log2_data_nsp)
        normalized_nsp_log2 = normalize(log2_data_nsp, medians_light)
        ic(normalized_nsp_log2)
        normalized_nsp = np.power(2, normalized_nsp_log2) # Inverse log2 transformation
        ic(normalized_nsp)
        plot_normalized_distributions(normalized_nsp, 'Normalized Distributions of Protein Expression Across NSP Samples')
    
    # Combine the normalized data from both datasets
    total_normalized = normalized_light.add(normalized_nsp, fill_value=0)
    
    # Saving the DataFrames to CSV files
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    ic(normalized_light)
    ic(normalized_nsp)
    ic(total_normalized)
    save_to_csv(normalized_light, os.path.join(save_dir, 'light.csv'))
    save_to_csv(normalized_nsp, os.path.join(save_dir, 'nsp.csv'))
    save_to_csv(total_normalized, os.path.join(save_dir, 'total.csv'))
    



if __name__ == '__main__':
    path = 'G:/My Drive/Data/data/testing pipeline dev/eif4f/'
    normalize_samples(path)













# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# import os

# def load_and_prepare_data(file_path):
#     """Load data from a CSV file and prepare it for analysis."""
#     try:
#         data = pd.read_csv(file_path)
#         data.replace([np.inf, -np.inf], np.nan, inplace=True)
#         data.dropna(how='all', subset=data.columns[1:], inplace=True)
#         return data
#     except Exception as e:
#         print(f"Error loading data: {e}")
#         return None

# def plot_log2_distributions(data, title):
#     """Plot log2 distributions of the data."""
#     log2_data = np.log2(data.iloc[:, 1:])
#     plt.figure(figsize=(20, 10))
#     sns.boxplot(data=log2_data)
#     plt.xticks(rotation=90)
#     plt.ylabel('Log2 Expression')
#     plt.title(title)
#     plt.tight_layout()
#     plt.show()
#     return log2_data

# def normalize_and_convert(log2_data, medians):
#     """Normalize the log2 data and convert back to original scale."""
#     normalized_log2_data = log2_data.subtract(medians)
#     normalized_data = np.power(2, normalized_log2_data) # Inverse log2 transformation
#     return normalized_data

# def plot_normalized_distributions(data, title):
#     """Plot distributions of the normalized data."""
#     plt.figure(figsize=(20, 10))
#     sns.boxplot(data=data)
#     plt.xticks(rotation=90)
#     plt.ylabel('Normalized Expression')
#     plt.title(title)
#     plt.tight_layout()
#     plt.show()

# def save_to_csv(df, file_path):
#     """Save DataFrame to a CSV file."""
#     try:
#         df.to_csv(file_path, index=False)
#         print(f"Saved to {file_path}")
#     except Exception as e:
#         print(f"Error saving file: {e}")

# # Directories and file paths
# base_dir = 'G:/My Drive/Data/data/eIF4F pilot/imputed'
# file_path_light = os.path.join(base_dir, 'light.csv')
# file_path_nsp = os.path.join(base_dir, 'nsp.csv')
# save_dir = os.path.join(base_dir, 'normalized')

# # Load and process the first dataset
# data_light = load_and_prepare_data(file_path_light)
# if data_light is not None:
#     log2_data_light = plot_log2_distributions(data_light, 'Log2 Distributions of Protein Expression Across Light Samples')
#     medians_light = log2_data_light.median()
#     normalized_light = normalize_and_convert(log2_data_light, medians_light)
#     plot_normalized_distributions(normalized_light, 'Normalized Distributions of Protein Expression Across Light Samples')

# # Load and process the second dataset
# data_nsp = load_and_prepare_data(file_path_nsp)
# if data_nsp is not None:
#     log2_data_nsp = plot_log2_distributions(data_nsp, 'Log2 Distributions of Protein Expression Across NSP Samples')
#     normalized_nsp = normalize_and_convert(log2_data_nsp, medians_light)
#     plot_normalized_distributions(normalized_nsp, 'Normalized Distributions of Protein Expression Across NSP Samples')

# # Combine the normalized data from both datasets
# total_normalized = normalized_light.add(normalized_nsp, fill_value=0)

# # Saving the DataFrames to CSV files
# if not os.path.exists(save_dir):
#     os.makedirs(save_dir)

# save_to_csv(normalized_light, os.path.join(save_dir, 'light.csv'))
# save_to_csv(normalized_nsp, os.path.join(save_dir, 'nsp.csv'))
# save_to_csv(total_normalized, os.path.join(save_dir, 'total.csv'))

