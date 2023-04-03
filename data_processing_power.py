#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 19:14:28 2023

@author: jha3
"""

# In[]
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
from datetime import date, datetime, timedelta, time

# In[]
from data_processing_base import *

# In[]
def get_power_data(nc_data):   
    power_ts = np.array(nc_data['POWERTS'])
    power_avg = np.array(nc_data['POWERTS'].mean(dim = 'Time'))
    
    return power_ts, power_avg


# In[]
def create_power_data (nc_data, start_time):
    pickled_data = create_base_data(nc_data, start_time)
    
    # Pickled power data
    pickled_data['power_inst'], pickled_data['power_avg'] = get_power_data(nc_data)
    
    return pickled_data


# In[]:
def combine_power_for_intervals_of_case(pickled_data_combined, pickled_data_for_interval, compute_power_stdev):
    power_avg_key = 'power_avg'
    power_stdev_key = 'power_stdev'
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
        
    # Append the stdev of power
    if compute_power_stdev:
        powert_stdev_for_interval = np.std(pickled_data_for_interval[power_inst_key])
        if power_stdev_key not in pickled_data_combined:
            pickled_data_combined[power_stdev_key] = powert_stdev_for_interval
        else:
            pickled_data_combined[power_stdev_key] = np.vstack((pickled_data_combined[power_stdev_key], powert_stdev_for_interval))
    
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