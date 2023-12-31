# -*- coding: utf-8 -*-
"""
Step 1: Create metadata

Adjust the path so that the metadata is created based on the Runs in the noralized intensities

then use this metadata in the next steps

"""

import pandas as pd
path = "G:/My Drive/Data/data/spikein data/protein intensities/"
df = pd.read_csv(f'{path}light_href.csv')
runs = df.columns.values.tolist()[1:]
print(runs)
# meta_dict = {'Run':runs,
# 'Sample':['FAC_1', 'FAC_2', 'FAC_3', 'Gu3340_1', 'Gu3340_2', 'Gu3340_3','both_1', 'both_2', 'both_3','control_1', 'control_2', 'control_3'],
#   'Treatment':['FAC','FAC','FAC', 'Gu3340', 'Gu3340', 'Gu3340', 'both', 'both', 'both', 'control', 'control', 'control']}

# meta_df = pd.DataFrame(meta_dict)
# meta_df.to_csv(f'{path}meta.csv')

# import pandas as pd
# path = "G:/My Drive/Data/data/eIF4F optimization/protein intensities/"
# df = pd.read_csv(f'{path}light_href.csv')
# runs = df.columns.values.tolist()[1:]
# meta_dict = {'Run':runs,
# 'Sample':['eIF4E+ 2h 1', 'eIF4E+ 2h 2', 'eIF4G1+ 2h 1', 'eIF4G1+ 2h 2', 'eIF4G2+ 2h 1', 'eIF4G2+ 2h 2','eIF4G3+ 2h 1', 'eIF4G3+ 2h 2',
#           'eIF4E- 2h 1', 'eIF4E- 2h 2', 'eIF4G1- 2h 1', 'eIF4G1- 2h 2', 'eIF4G2- 2h 1', 'eIF4G2- 2h 2','eIF4G3- 2h 1', 'eIF4G3- 2h 2',
#           'eIF4E+ 4h 1', 'eIF4E+ 4h 2', 'eIF4G1+ 4h 1', 'eIF4G1+ 4h 2', 'eIF4G2+ 4h 1', 'eIF4G2+ 4h 2','eIF4G3+ 4h 1', 'eIF4G3+ 4h 2',
#           'eIF4E- 4h 1', 'eIF4E- 4h 2', 'eIF4G1- 4h 1', 'eIF4G1- 4h 2', 'eIF4G2- 4h 1', 'eIF4G2- 4h 2','eIF4G3- 4h 1', 'eIF4G3- 4h 2',
#           'eIF4E+ 8h 1', 'eIF4E+ 8h 2', 'eIF4G1+ 8h 1', 'eIF4G1+ 8h 2', 'eIF4G2+ 8h 1', 'eIF4G2+ 8h 2','eIF4G3+ 8h 1', 'eIF4G3+ 8h 2',
#           'eIF4E- 8h 1', 'eIF4E- 8h 2', 'eIF4G1- 8h 1', 'eIF4G1- 8h 2', 'eIF4G2- 8h 1', 'eIF4G2- 8h 2','eIF4G3- 8h 1', 'eIF4G3- 8h 2']}

# meta_df = pd.DataFrame(meta_dict)
# meta_df.to_csv(f'{path}meta.csv', index=False)