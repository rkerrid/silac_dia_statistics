# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:15:50 2023

@author: robbi
"""

# get control and treatment pairs for light data based on metadata
# where >=2 ids per pg in both pairs, make a distribution
# shift = median treated - median control
# normalizing the treated to the control is treated + shift for all non NaN values in treated
# output normalized light df + shift values to be applied to NSPdf

# import NSP and shift vals
# for all NSP treated, add the shift value for all non NaN values
# output normalized NSP data

# import light and NSP unormalized, plot histograms of all non NaN vals, colour sample/treatmentwise
# do the same for the normalized and compare

