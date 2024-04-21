# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:48:39 2024

@author: robbi
"""

import alphastats 
import pandas as pd
import matplotlib.pyplot as plt  # Make sure matplotlib.pyplot is imported as plt
import plotly.express as px


def generate_alphastats_objects(path, meta):
    meta_file = pd.read_csv(meta, sep=',')
    path = f'{path}/protein_groups/'
    intensity_cols = meta_file['Sample'].values.tolist()
    # loader_total = alphastats.GenericLoader(f"{path}total.csv", 
    #                                   intensity_column = intensity_cols,
    #                                     index_column="Protein.Group")

    loader_nsp = alphastats.GenericLoader(f"{path}nsp.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

    loader_light = alphastats.GenericLoader(f"{path}light.csv", 
                                      intensity_column = intensity_cols,
                                        index_column="Protein.Group")

    # df_total = alphastats.DataSet(
    #     loader = loader_total,
    #     metadata_path = meta,
    #     sample_column = 'Sample')

    df_light = alphastats.DataSet(
        loader = loader_light,
        metadata_path = meta,
        sample_column = 'Sample')

    df_nsp = alphastats.DataSet(
        loader = loader_nsp,
        metadata_path = meta,
        sample_column = 'Sample')


    
    return df_nsp, df_light  #, df_total

meta = 'G:/My Drive/Data/data/20240410 SRP AID/meta.csv'
path = 'G:/My Drive/Data/data/20240410 SRP AID/'

nsp, light = generate_alphastats_objects(path, meta)
light.preprocess()
df = light.rawmat
# filtered_df = df[df['column_name'].str.contains('P61011')]
print(df.columns)


filtered_columns = [col for col in df.columns if 'P61011' in col]

# Filter the DataFrame to keep only these columns
filtered_df = df[filtered_columns]
# Now transpose the DataFrame
transposed_df = filtered_df.transpose()

# Create a bar plot
plt.figure(figsize=(10, 8))  # Adjust size as needed
bar_plot = transposed_df.plot(kind='bar', legend=False, ax=plt.gca())  # Use gca to get the current axes for plt

transposed_df.plot(kind='bar', legend=False)  # You can turn legend on by setting legend=True
plt.title('Protein Group P61011 Concentrations')
plt.ylabel('Log2-transformed Concentration')
plt.xlabel('Samples/Measurements')
plt.xticks(rotation=45)  # Rotates labels to avoid overlap
bar_plot.set_xticklabels(transposed_df.index, rotation=45, ha="right")  # Adjust rotation and alignment for better visibility

plt.tight_layout()  # Adjust layout to make room for label rotation

# Show plot
plt.show()
