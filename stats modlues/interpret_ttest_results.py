# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:07:03 2023

@author: robbi
"""

import pandas as pd
import numpy as np
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic
import plotly.express as px
import ipywidgets as widgets
from IPython.display import display
import os
from manage_directories import create_directory

global significant_proteins_df
significant_proteins_df = None
sig = None
fc = None

def scan_folder(path):
    # List to hold file names
    file_names = []

    # Check if the path exists and is a directory
    if os.path.exists(path) and os.path.isdir(path):
        # Iterate over the files in the directory
        for entry in os.listdir(path):
            full_path = os.path.join(path, entry)
            # Check if it's a file and not a directory
            if os.path.isfile(full_path):
                file_names.append(entry)
    else:
        print("The provided path does not exist or is not a directory.")

    return file_names


def create_volcano_plot(df, title, pois, path):
    def volcano(significance_level=0.05, fold_change=1.0, show_labels=False):
        df['-log10(p_value)'] = -np.log10(df['pval'])

        # Determine the color and labels for each point
        df['color'] = 'blue'  # Default color
        significant = df['-log10(p_value)'] >= -np.log10(significance_level)
        high_fold_change = np.abs(df['log2fc']) >= fold_change
        df['is_significant'] = significant & high_fold_change
        df.loc[df['is_significant'], 'color'] = 'red'
        
        df['protein_name'] = df['Protein.Group'].apply(lambda x: x.split('-')[1] if '-' in x else x)


        if show_labels:
            df['labels'] = np.where(df['is_significant'], df['protein_name'], None)
        else:
            df['labels'] = None
        if len(pois) >0:
            df['labels'] = np.where(df['protein_name'].isin(pois), df['protein_name'], None)
            
        fig = px.scatter(df, x='log2fc', y='-log10(p_value)', hover_name='protein_name', color='color',
                         title=title,
                         labels={'log2fc': 'Log2 Fold Change', '-log10(p_value)': '-log10(p-value)'},
                         color_discrete_map={'red': 'red', 'blue': 'blue'},
                         text='labels')

        significant_threshold = -np.log10(significance_level)
        fig.add_hline(y=significant_threshold, line_dash="dash", line_color="red")

        # Vertical lines for fold change
        fig.add_vline(x=fold_change, line_dash="dash", line_color="blue")
        fig.add_vline(x=-fold_change, line_dash="dash", line_color="blue")

        fig.update_layout(height=600, width=800)
        fig.update_traces(textposition='top center', textfont={'color':'black', 'size':12, 'family':'Arial, sans-serif'})
        sig = significance_level
        fc = fold_change
        return fig

    # def generate_significant_proteins_list():
    #     significant = df['-log10(p_value)'] >= -np.log10(sig_slider.value)
    #     high_fold_change = np.abs(df['log2fc']) >= fc_slider.value
    #     significant_proteins = df[significant & high_fold_change]

    #     significant_proteins['upregulated'] = significant_proteins['log2fc'] > 0
    #     significant_proteins['downregulated'] = significant_proteins['log2fc'] < 0

    #     result_df = significant_proteins[['Protein.Group', 'upregulated', 'downregulated']]
    #     result_df.columns = ['Proteins', 'Upregulated', 'Downregulated']

    #     return result_df
    def generate_significant_proteins_list():
        significant = df['-log10(p_value)'] >= -np.log10(sig_slider.value)
        high_fold_change = np.abs(df['log2fc']) >= fc_slider.value
        significant_proteins = df[significant & high_fold_change]
    
        # Create a new column for regulation status
        significant_proteins['Regulation'] = np.where(significant_proteins['log2fc'] > 0, 'Upregulated', 'Downregulated')
    
        # Create the result DataFrame with two columns
        result_df = significant_proteins[['Protein.Group', 'Regulation']]
        result_df['protein_name'] = result_df['Protein.Group'].apply(lambda x: x.split('-')[1] if '-' in x else x)
        # result_df.columns = ['Proteins', 'Regulation']
        
        # Sort by Regulation
        result_df = result_df.sort_values(by='Regulation')
        # Save to CSV
        create_directory(path, 'metascape_lists')
        # create_directory(f'{path}metascape_lists', '{title}_l2fc:{fc}_sig:{sig}')
        result_df.to_csv(f'{path}metascape_lists/{title}.csv', index=False)
    
        return result_df

# Example usage
# df = your_dataframe
# sig_slider_value = value_from_slider
# fc_slider_value = value_from_slider
# output_filename = "significant_proteins.csv"
# generate_significant_proteins_list(df, sig_slider_value, fc_slider_value, output_filename)


    def update_plot(change):
        with output:
            output.clear_output(wait=True)
            display(volcano(sig_slider.value, fc_slider.value, label_button.value))

    def on_button_clicked(b):
        global significant_proteins_df
        significant_proteins_df = generate_significant_proteins_list()
        with output_df:
            output_df.clear_output(wait=True)
            display(significant_proteins_df)

    # Sliders and Buttons
    sig_slider = widgets.FloatSlider(value=0.05, min=0.001, max=0.1, step=0.01, description='Significance Level:', continuous_update=False)
    fc_slider = widgets.FloatSlider(value=1.0, min=0.1, max=5.0, step=0.1, description='Fold Change:', continuous_update=False)
    label_button = widgets.ToggleButton(value=False, description='Show Labels', icon='eye')
    generate_button = widgets.Button(description="Generate List")

    # Observers
    sig_slider.observe(update_plot, names='value')
    fc_slider.observe(update_plot, names='value')
    label_button.observe(update_plot, names='value')
    generate_button.on_click(on_button_clicked)

    # Output Widgets
    output = widgets.Output()
    output_df = widgets.Output()

    # Display Widgets
    with output:
        display(volcano(sig_slider.value, fc_slider.value, label_button.value))

    control_widgets = widgets.VBox([sig_slider, fc_slider, label_button, generate_button])
    display(control_widgets, output, output_df)


def reduce_df_size(df, top_n):
    
    # Sort the DataFrame by the p-value column in ascending order and select the top x rows
    filtered_df = df.sort_values(by='pval', ascending=True).head(top_n)
    return filtered_df

def loop_and_plot_results(path, result_list, pois):
    for result in result_list:
        name = result[:-4]
        df = pd.read_csv(f'{path}/{result}')
        if len(df) > 1000:
            df = reduce_df_size(df, 1000)
        create_volcano_plot(df, name, pois, path)
    
    
    
    
    
    
    
    