
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

# pd.set_option('display.max_rows', None)
ic.disable()

def import_metadata(path):
    """
    Import metadata from a given path.
    """
    return pd.read_csv(f"{path}")

def subset_metadata(metadata, subset):
    """
    Subset the metadata based on the given subset.
    """
    if subset:
        metadata = metadata[metadata['Treatment'].isin(subset)]
    return metadata

def load_dataframes(path):
    """
    Load dataframes from given path and quantification method.
    """
    path = f'{path}/normalized/'
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

def add_lowest_sample_mean_for_row(df, metadata): # should this function just take lowest mean for control? Or median of at least two valid values per sample?
    df['lowest_mean'] = np.nan

    for group in metadata['Treatment'].unique():
        sub_meta = metadata[metadata['Treatment'] == group]
        cols = sub_meta['Sample'].tolist()

        for index, row in df.iterrows():
            valid_values = row[cols].dropna()
            if len(valid_values) == len(cols): # if there are at least two value to obtain a mean for this protein
                row_mean = valid_values.mean()
                if pd.isna(df.at[index, 'lowest_mean']) or row_mean < df.at[index, 'lowest_mean']:
                    df.at[index, 'lowest_mean'] = row_mean

    return df

def get_global_imputation_values(df, metadata, channel, control_samples):
    df_copy = df.copy(deep = True)
    df_copy = subset_data(df_copy, metadata)
    cols = df_copy.columns[1:]
    df_copy = df_copy[cols]
    # just get the control nsp distribution if the data is the nsp dataset
    subset = ""
    if channel == 'nsp':
        cols = control_samples
        df_copy = df_copy[cols] # need to make dynamic
        subset = 'control samples'
        
    else:
        subset = 'all samples'
    ic(subset)
    data = df_copy.dropna()
    ic(data)
    mu, std = stats.norm.fit(data)
    ic(mu)
    ic(std)

    global_mu = mu - 1.8 * std
    global_std = std * 0.3
    ic(global_mu)
    ic(global_std)
    
    # plot distribution of global dataset
    # Creating the histogram
    # Flatten the DataFrame to a single series
    flat_data = data[cols].values.flatten()
    plt.hist(flat_data, bins=300)
    plt.axvline(global_mu, color='green', linestyle='dashed', linewidth=2)
    plt.axvline(mu, color='red', linestyle='dashed', linewidth=2)


    # Adding titles and labels
    plt.title(f'Histogram of {channel} {subset}')
    plt.xlabel('Frequency')
    plt.ylabel('Log2 intensity')
    plt.show()
    return global_mu, global_std



def impute_values(df, metadata, channel):
    """
    Impute missing values in the dataframe based on metadata.
    """
    df_copy = df.copy(deep=True)
    df_copy = subset_data(df_copy, metadata)
    cols = df_copy.columns.values.tolist()

    global_mu = df['global_mu'].iloc[0]
    global_std = df['global_std'].iloc[0]

    # Create a DataFrame to keep track of NaNs before imputation
    nan_before_imputation = df[cols].isna()

    # Vectorized calculation for adjusted mean
    adjusted_mu = np.where(df['lowest_mean'] + global_std < global_mu, 
                           df['lowest_mean'] - global_std, 
                           global_mu)

    # Vectorized imputation for each column
    for col in cols:
        # Identify rows with NaN values for this column
        nan_rows = df[col].isna()

        # Generate random values only for the NaN values
        random_values = np.random.normal(adjusted_mu[nan_rows], global_std, size=nan_rows.sum())

        # Update only NaN values in the original DataFrame
        df.loc[nan_rows, col] = random_values

    # Update the 'mu_used_for_imputation' column
    df['mu_used_for_imputation'] = adjusted_mu

    # Update 'was_imputed' to True where NaNs were present before imputation
    df['was_imputed'] = nan_before_imputation.any(axis=1)
    
    print('Inspect df after all annotation')
    df.to_csv('G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/H/protein_groups_statistics/normalized/imputed/annotated.csv', sep=',')
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
    
    # ic(imputed_vals.columns[1:])
    for col in imputed_vals.columns[1:]:
        plt.figure(figsize=(10, 6))

        sns.histplot(imputed_vals[col], color='blue', label='Imputed', alpha=0.6, bins=300)
        sns.histplot(non_imputed[col], color='orange', label='Original', alpha=0.6, bins=300)
        plt.title(f"{title} - {col}")
        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.legend()
        plt.show()

def annotate_df(df, metadata, channel, control_samples):
    print('annotate df with the lowest mean value of observed proteins within each sample group per row')
    global_mu, global_std = get_global_imputation_values(df, metadata, channel, control_samples)
    df['global_mu'] = global_mu
    df['global_std'] = global_std
    df['mu_used_for_imputation'] = np.nan
        
    df = add_lowest_sample_mean_for_row(df, metadata)
        
    return df
    
def process_intensities(path, control_samples, meta, subset=[], plot_imputation=False):
    """
    Main function to process protein intensities.
    """
    print('import data')
    path = f'{path}/statistics/'
    metadata = import_metadata(meta)
    light, nsp = load_dataframes(path)
    print('Subset data')
    metadata = subset_metadata(metadata, subset)
    light = subset_data(light, metadata)
    nsp = subset_data(nsp, metadata)
    print('preprocess data')
    light = preprocess_dataframe(light)
    nsp = preprocess_dataframe(nsp)
    print('annotate df with std,mu,lowestval, etc')
    # for nsp data I need the control samples from the control_samples list
    # Filtering and obtaining the samples
    matching_control_samples = metadata[metadata['Treatment'].isin(control_samples)]['Sample']
    
    # Convert to a list
    control_samples = matching_control_samples.tolist()
    ic(control_samples)
    
    light_annotated = annotate_df(light, metadata, 'light', control_samples)
    nsp_annotated = annotate_df(nsp, metadata, 'nsp', control_samples)
    light_annotated_copy = light_annotated.copy(deep = True)
    nsp_annotated_copy = nsp_annotated.copy(deep = True)
    print('preform imputation')
    light = impute_values(light_annotated, metadata,'light')
    nsp = impute_values(nsp_annotated, metadata, 'nsp')

    if plot_imputation:
        plot_histograms(light, 'Light', metadata)
        plot_histograms(nsp, 'Nsp', metadata)
        
    light = subset_data(light, metadata)
    nsp = subset_data(nsp, metadata)
    
    #base 2 exponentiation before saving
    light.iloc[:,1:] = 2**light.iloc[:,1:]
    nsp.iloc[:,1:] = 2**nsp.iloc[:,1:]
    
    save_imputed_data(light, nsp, metadata, path)
    return light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy

def save_imputed_data(light_df, nsp_df, metadata, path):
    """
    Save the imputed data to a specified path.
    """
    create_directory(f"{path}", "imputed")
    light_df.to_csv(f"{path}imputed/light.csv", sep=',', index=False)
    nsp_df.to_csv(f"{path}imputed/nsp.csv", sep=',', index=False)
    metadata.to_csv(f"{path}imputed/meta.csv", sep=',', index=False)
    
# if __name__ == '__main__':
#     path = 'G:/My Drive/Data/data/poc4/H/protien_groups/normalized/'
#     light, nsp, light_annotated, nsp_annotated, light_annotated_copy, nsp_annotated_copy = process_intensities(path, plot_imputation=True)

