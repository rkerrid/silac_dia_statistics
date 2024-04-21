# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 13:47:56 2024

@author: robbi
"""

from silac_dia_tools.pipeline_refactored.pipeline import Pipeline as pileline
from sdia_stats.preprocessing import adapted_imputation 
from sdia_stats.preprocessing import normalize_samples
from sdia_stats.preprocessing import filter_contams_and_non_valid_rows
from sdia_stats.qc import normalization_qc
from sdia_stats.statistics import ttest
from sdia_stats.visualization.interpret_ttest_results import loop_and_plot_results
import pandas as pd

meta = 'G:/My Drive/Data/data/20240410 SRP AID/meta.csv'
path = 'G:/My Drive/Data/data/20240410 SRP AID/'

normalization_qc.pca(path,meta)