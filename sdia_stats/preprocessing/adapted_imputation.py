
"""
Created on Fri Jan 12 11:58:31 2024
@author: robbi
"""
from tqdm import tqdm
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from sdia_stats.utils.manage_directories import create_directory
from icecream import ic
ic.enable()

def import_metadata(path):
    """
    Import metadata from a given path.
    """
    return pd.read_csv(f"{path}meta.csv")

def subset_metadata(metadata, subset):
    """
    Subset the metadata based on the given subset.
    """
    if subset:
        metadata = metadata[metadata['Treatment'].isin(subset)]
    return metadata

def load_dataframes(path, quantification):
    """
    Load dataframes from given path and quantification method.
    """
    light = pd.read_csv(f"{path}light.csv", sep=',')
    nsp = pd.read_csv(f"{path}nsp.csv", sep=',')
    return light, nsp

def preprocess_dataframe(df):
    """
    Preprocess the dataframe by replacing infinite values and thresholding.
    """
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.iloc[:, 1:] = df.iloc[:, 1:].apply(lambda x: np.where(x < 0.001, np.nan, x))
    # log2 data
    df.iloc[:, 1:] = np.log2(df.iloc[:, 1:].replace([np.inf, -np.inf], np.nan))

    return df

def subset_data(df, metadata):
    """
    Subset the dataframe based on the metadata.
    """
    relevant_samples = metadata['Sample'].values.tolist()
    columns_to_keep = ['Protein.Group'] + relevant_samples
    return df[columns_to_keep]

def filter_valid_values(df, metadata):
    """
    Filter the dataframe for valid values based on metadata.
    """
    df['keep_row'] = df.apply(lambda row: row[metadata['Sample'].tolist()].notna().sum() >= 2, axis=1)
    return df[df['keep_row']].drop('keep_row', axis=1)

def add_lowest_sample_mean_for_row(df, metadata):
    df['lowest_mean'] = np.nan

    for group in metadata['Treatment'].unique():
        sub_meta = metadata[metadata['Treatment'] == group]
        cols = sub_meta['Sample'].tolist()

        for index, row in df.iterrows():
            valid_values = row[cols].dropna()
            if len(valid_values) == len(cols):
                row_mean = valid_values.mean()
                if pd.isna(df.at[index, 'lowest_mean']) or row_mean < df.at[index, 'lowest_mean']:
                    df.at[index, 'lowest_mean'] = row_mean

    return df

def get_global_imputation_values(df, metadata):
    df_copy = df.copy(deep = True)
    df_copy = subset_data(df_copy, metadata)
    cols = df_copy.columns[1:]
    df_copy = df_copy[cols]
    data = df_copy.dropna()
    ic(data)
    global_mu, global_std = stats.norm.fit(data)
    ic(global_mu)
    ic(global_std)

    global_mu = global_mu - 1.8 * global_std
    global_std = global_std * 0.3
    ic(global_mu)
    ic(global_std)
    
    return global_mu, global_std

def impute_values(df, metadata):
    """
    Impute missing values in the dataframe based on metadata.
    """
    
    df_copy = df.copy(deep=True)
    df_copy = subset_data(df_copy, metadata)
    cols = df_copy.columns.values.tolist()
    
    df['was_imputed'] = False
    global_mu = df['global_mu'].iloc[0]
    global_std = df['global_std'].iloc[0]
    
    
    for col in tqdm(cols, desc="Imputing columns"):
        for index, row in df.iterrows():
            if pd.isna(row[col]):
                if row['lowest_mean'] + global_std*0.5 < global_mu:
                    adjusted_mu = row['lowest_mean'] - global_std
                    imputed_value = np.random.normal(adjusted_mu, global_std)
                else:
                    imputed_value = np.random.normal(global_mu, global_std)
                df.at[index, 'mu_used_for_imputation'] = imputed_value
                df.at[index, col] = imputed_value
                df.at[index, 'was_imputed'] = True

    return df

def plot_histograms(df, title, metadata):
    """
    Plot histograms for the dataframe.
    """
    df = df.copy(deep=True)
    imputed_vals = df[df['was_imputed'] == True]
    non_imputed = df[df['was_imputed'] == False]
    imputed_vals = subset_data(imputed_vals, metadata)
    non_imputed = subset_data(non_imputed, metadata)
    
    ic(imputed_vals.columns[1:])
    for col in imputed_vals.columns[1:]:
        plt.figure(figsize=(10, 6))

        sns.histplot(imputed_vals[col], color='blue', label='Imputed', alpha=0.6, bins=300)
        sns.histplot(non_imputed[col], color='orange', label='Original', alpha=0.6, bins=300)
        plt.title(f"{title} - {col}")
        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.legend()
        plt.show()

def annotate_df(df, metadata):
    print('annotate df with the lowest mean value of observed proteins within each sample group per row')
    global_mu, global_std = get_global_imputation_values(df, metadata)
    df['global_mu'] = global_mu
    df['global_std'] = global_std
    df['mu_used_for_imputation'] = np.nan
    df = add_lowest_sample_mean_for_row(df, metadata)
    return df
    
def process_intensities(path, subset=[], plot_imputation=False, quantification='href'):
    """
    Main function to process protein intensities.
    """
    print('import data')
    metadata = import_metadata(path)
    light, nsp = load_dataframes(path, quantification)
    print('Subset data')
    metadata = subset_metadata(metadata, subset)
    light = subset_data(light, metadata)
    nsp = subset_data(nsp, metadata)
    print('preprocess data')
    light = preprocess_dataframe(light)
    nsp = preprocess_dataframe(nsp)
    print('filter for valid values')
    light = filter_valid_values(light, metadata)
    nsp = filter_valid_values(nsp, metadata)
    print('annotate df with std,mu,lowestval, etc')
    light_annotated = annotate_df(light, metadata)
    nsp_annotated = annotate_df(nsp, metadata)
    print('preform imputation')
    light = impute_values(light_annotated, metadata)
    nsp = impute_values(nsp_annotated, metadata)

    if plot_imputation:
        plot_histograms(light, 'Light', metadata)
        plot_histograms(nsp, 'Nsp', metadata)
        
    light = subset_data(light, metadata)
    nsp = subset_data(nsp, metadata)
    
    #base 2 exponentiation before saving
    light.iloc[:,1:] = 2**light.iloc[:,1:]
    nsp.iloc[:,1:] = 2**nsp.iloc[:,1:]
    
    save_imputed_data(light, nsp, path)
    return light, nsp, light_annotated, nsp_annotated

def save_imputed_data(light_df, nsp_df, path):
    """
    Save the imputed data to a specified path.
    """
    create_directory(f"{path}", "imputed")
    light_df.to_csv(f"{path}imputed/light.csv", sep=',', index=False)
    nsp_df.to_csv(f"{path}imputed/nsp.csv", sep=',', index=False)

if __name__ == '__main__':
    path = 'G:/My Drive/Data/data/poc4/H/normalized/'
    light, nsp, light_annotated, nsp_annotated = process_intensities(path, plot_imputation=True)


# from sdia_stats.utils.manage_directories import create_directory

# import pandas as pd
# import numpy as np
# from scipy import stats
# import seaborn as sns
# import matplotlib.pyplot as plt

# # from icecream import ic

#  # import metadata
 

# from icecream import ic
# ic.disable()
#  # import metadata
 

# def import_meta(path):
#      metadata = pd.read_csv(f"{path}meta.csv")
#      return metadata
 
# def get_dataframes(path, quantification):
#     # total = pd.read_csv(f"{path}protein intensities/total_{quantification}.csv", sep=',')
#     light = pd.read_csv(f"{path}light.csv", sep=',')
#     nsp = pd.read_csv(f"{path}nsp.csv", sep=',')
#     return  light, nsp
     
# def replace_values(df):
#     df.replace([np.inf, -np.inf], np.nan, inplace=True)
#     df.iloc[:,1:] = df.iloc[:,1:].apply(lambda x: np.where(x < 0.001, np.nan, x))
#     return df

         
# def filter_for_valid_values(df, metadata): 
#      df['keep_row'] = False
   
#      for group in metadata['Treatment'].unique():
#          sub_meta = metadata[metadata['Treatment'] == group]
#          cols = sub_meta['Sample'].tolist()
#          # Check that we have at least 2 columns to compare, if not continue to next group
#          if len(cols) < 2:
#              continue
 
#          # Calculate the sum of valid (not NaN) values for the columns in cols
#          valid_sum = df[cols].notna().sum(axis=1)

#          # Update 'keep_row' to True where at least 2 of the columns are not NaN
#          df.loc[valid_sum >= 2, 'keep_row'] = True
#      df = df[df['keep_row']].copy()

#      df = df[df['keep_row']]
#      df.drop('keep_row', axis=1, inplace=True)
#      return df     
 
# # imputation
# def create_distribution(data):
#     mu, std = stats.norm.fit(data)
#     mu = mu - 1.8 * std
#     std = 0.25*std
   
#     return mu, std
    
# def impute(x, mu, std): # at this point x is a row[col]
    
#     # if imputed_only:
#     if np.isnan(x):
#         x = np.random.normal(mu, std)
#         return x, True
#     else:
#         return x, False
        



# def add_lowest_sample_mean_for_row(df, metadata):
#     # Extract the sample columns (assuming they start with 'sample_')
#     df['lowest_mean'] = np.nan
#     ic(df)
#     for group in metadata['Treatment'].unique():
#         sub_meta = metadata[metadata['Treatment'] == group]
#         cols = sub_meta['Sample'].tolist()
#         # now itterate through sub df and if there are at least 2 valid vals, return the mean
#         # if the current mean is not NaN or less than the new mean, the new lowest mean is the mean of that sample group for that protein
#         for index, row in df.iterrows():
#             new_mean = 0
#             # Calculate the sum of valid (not NaN) values for the columns in cols
#             valid_sum = row[cols].notna().sum()
#             ic(valid_sum)
#             if valid_sum== len(cols):
#                 new_mean = row[cols].mean()
#                 ic(new_mean)
#             else:
#                 new_mean = 10000000000000000000000000000000000000000000000000000000000000 # fix this
#             if not np.isnan(row['lowest_mean']) or row['lowest_mean'] > new_mean:
#                     ic(row['lowest_mean'])
#                     row['lowest_mean'] = new_mean
#                     ic(row['lowest_mean'])
#     return df




# def perform_imputation(df, metadata):
#     cols = df.columns.values.tolist()[1:]
#     df['was_imputed'] = False
#     df = add_lowest_sample_mean_for_row(df, metadata)
    
#     for condition in cols:
#         data = df[condition].dropna()
#         global_mu, global_std = create_distribution(data)
        
#         for index, row in df.iterrows():
#             if row['lowest_mean'] < global_mu + global_std:
#                 adjusted_mu = row['lowest_mean'] - global_std
#                 imputed_value, was_imputed = impute(row[condition], adjusted_mu, global_std)
#             else:
#                 imputed_value, was_imputed = impute(row[condition], global_mu, global_std)
            
#             df.at[index, condition] = imputed_value
#             if was_imputed and df.at[index, 'was_imputed'] == True:
#                 df.at[index, 'was_imputed'] = True

   
    
#     return df






# # def plot_histogram(df, imputed_values, title):
# #   # Plot histograms for each column in df1 and df2 on the same plot with different colors
# #     # Create a single figure with 2 rows and 3 columns
  
# #     # Plot histograms for each column in df and imputed_values
# #     for  col in df.columns.values.tolist()[1:]:
# #         imputed_data = imputed_values[col].dropna()
# #         original_data = df[col].dropna()
        
# #         plt.hist(original_data, bins=200, alpha=0.5, label='original data', color='blue')
# #         plt.hist(imputed_data, bins=200, alpha=0.5, label='imputed values', color='green')
# #         plt.title(f'{title} {col} Histogram')
# #         plt.xlabel('Log2 intensity')
# #         plt.ylabel('Frequency')
# #         # plt.legend()

  

# #         # Display the figure
# #         plt.show()

# def subset_data(df,  metadata):
#     # need to map the cols to keep witht he metadata and drop/save cols that meet re requirements

#     relevant_samples = metadata['Sample'].values.tolist()
 
#     columns_to_keep = ['Protein.Group'] + relevant_samples

#     df = df[columns_to_keep]
    
#     return df

# def subset_metadata(metadata, subset):
#     if len(subset)>0:
#         metadata = metadata[metadata['Treatment'].isin(subset)]
#     return metadata


# def plot_imputed_results(df, meta, title):
    
#     imputed_vals = df[df['was_imputed'] == True]
#     non_imputed = df[df['was_imputed'] == False]
#     for col in df.columns.values.tolist()[1:]:
#         sns.histplot(imputed_vals[col])
#         sns.histplot(non_imputed[col])
#         plt.show()


# def process_intensities(path, subset = [], plot_imputation=False, quantification='href'):
#     metadata = import_meta(path)
#     # groups = metadata[metadata_sample_group].unique()
    
#     light, nsp = get_dataframes(path, quantification)
#     print('imported')
#     metadata = import_meta(path)
#     metadata = subset_metadata(metadata, subset)
#     light = subset_data(light, metadata)
#     nsp = subset_data(nsp, metadata)
#     print('subseted')
#     # replace NaN and inf values
#     light = replace_values(light)
#     nsp = replace_values(nsp)
  
#     print('replace')
  
    
#     # filter for valid values
#     light = filter_for_valid_values(light, metadata)
#     nsp = filter_for_valid_values(nsp, metadata)
    
#     print('valid values')
  
#     # log transform
#     light.iloc[:,1:] = np.log2(light.iloc[:,1:])
#     nsp.iloc[:,1:] = np.log2(nsp.iloc[:,1:])
    
    
#     # impute with gausian shift
#     nsp_df = perform_imputation(nsp, metadata)
#     light_df = perform_imputation(light, metadata)
 
    
#     if plot_imputation:
#         plot_imputed_results(light_df,'Light')
#         plot_imputed_results(nsp_df, 'Nsp')
    
    
#     # base 2 exponentiation before saving
#     nsp_df.iloc[:,1:] = 2**nsp_df.iloc[:,1:]
#     light_df.iloc[:,1:] = 2**light_df.iloc[:,1:]
    
#     create_directory(f"{path}", "imputed")
#     metadata.to_csv(f"{path}/imputed/meta.csv", sep=',',index=False)
#     nsp_df.to_csv(f"{path}imputed/nsp.csv", sep=',',index=False)
#     light_df.to_csv(f"{path}imputed/light.csv", sep=',',index=False)
#     return light_df, nsp_df

# if __name__ == '__main__':
#     path = 'G:/My Drive/Data/data/poc4/H/protein intensities/'
#     path = 'G:/My Drive/Data/data/poc4/H/normalized/'
#     process_intensities(path, plot_imputation=True)
