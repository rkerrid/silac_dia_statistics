# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:50:09 2023

@author: robbi
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def test_function():
    print("hello world")
    
def plot_proteins(df, proteins, sample_columns):
    # Filter the DataFrame to include only the selected sample columns and the proteins of interest
    selected_df = df[["Protein.Group"] + sample_columns]
    selected_df = selected_df[selected_df["Protein.Group"].isin(proteins)]

    # Set the "Proteins" column as the index for easier plotting
    selected_df = selected_df.set_index("Protein.Group")

    # Create a plot with lines connecting data points for each sample
    plt.figure(figsize=(10, 6))
    for sample in sample_columns:
        plt.plot(selected_df.index, selected_df[sample], label=sample)

    # Customize the plot
    plt.xlabel("Proteins")
    plt.ylabel("Intensity")
    plt.title("Protein Intensity for Selected Proteins")
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    plt.legend()  # Add a legend to differentiate sample lines

    # Show the plot
    plt.show()

import protein_levels_plots as protplot
import pandas as pd
path = "G:/My Drive/Data/data/eIF4F optimization/protein intensities/"

df = pd.read_csv(f"{path}total_href.csv")
df['Protein.Group'] = df['Protein.Group'].str.split('-').str[1]
# df = df[df['Protein.Group']=="EIF4A1"]
# df = df[df['Protein.Group']=="EEFG1"]
# df = df[df['Protein.Group']=="EIF3A"]
# df = df[df['Protein.Group']=="PRAS40"]
# df = df[df['Protein.Group']=="EEF1D"]
# df = df[df['Protein.Group']=="EIF4A3"]
# df = df[df['Protein.Group']=="EIF4G1"]

# df = df[df['Protein.Group']=="EIF3B"]
# df = df[df['Protein.Group']=="EIF3C"]