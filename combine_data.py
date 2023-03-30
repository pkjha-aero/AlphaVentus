#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jha3
"""
# In[]
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
from datetime import date, datetime, timedelta, time

# In[]:
def combine_power_for_intervals_of_case(pickled_data_combined, pickled_data_for_interval):
    power_avg_key = 'power_avg'
    power_inst_key = 'power_inst'
    
    # Append the averaged power
    if power_avg_key not in pickled_data_combined:
        pickled_data_combined[power_avg_key] = pickled_data_for_interval[power_avg_key]
    else:
        pickled_data_combined[power_avg_key] = np.vstack((pickled_data_combined[power_avg_key], pickled_data_for_interval[power_avg_key]))
    
    # Append the instantaneous power
    if power_inst_key not in pickled_data_combined:
        pickled_data_combined[power_inst_key] = pickled_data_for_interval[power_inst_key]
    else:
        pickled_data_combined[power_inst_key] = np.vstack((pickled_data_combined[power_inst_key], pickled_data_for_interval[power_inst_key]))
    
    # Return the combined data
    return pickled_data_combined

# In[]:
def compute_power_pdf_combined (combined_power_data, num_bins = 20):
    
    hist, bin_edges = np.histogram(combined_power_data['power_inst']/1.0e6, bins = num_bins, density = True)
    bin_centers = bin_edges[:-1] + 0.5*np.diff(bin_edges)
    
    combined_power_data['bin_edges'] = bin_edges
    combined_power_data['bin_centers'] = bin_centers
    combined_power_data['hist'] = hist
    combined_power_data['pdf_integral'] = np.sum(hist * np.diff(bin_edges))
     
    return combined_power_data