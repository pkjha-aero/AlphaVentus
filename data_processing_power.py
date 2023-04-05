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
from scipy import signal
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


# In[]
def scale_power (pickled_power_data_interval, surplus_threshold = 4.0):
    power_avg_key = 'power_avg'
    power_inst_key = 'power_inst'
    pickled_power_data_interval_scaled = {}
    
    surplus = pickled_power_data_interval[power_inst_key]/1.0e6 - surplus_threshold
    surplus[np.where(surplus < 0)] = 0.0
    surplus_normalized = surplus/surplus.max()
    
    reduction = surplus*surplus_normalized
    
    pickled_power_data_interval_scaled[power_inst_key] = pickled_power_data_interval[power_inst_key] - reduction*1.0e6
    pickled_power_data_interval_scaled[power_avg_key] = pickled_power_data_interval_scaled[power_inst_key].mean()
    
    return pickled_power_data_interval_scaled

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
        power_stdev_for_interval = np.std(pickled_data_for_interval[power_inst_key])
        if power_stdev_key not in pickled_data_combined:
            pickled_data_combined[power_stdev_key] = power_stdev_for_interval
        else:
            pickled_data_combined[power_stdev_key] = np.vstack((pickled_data_combined[power_stdev_key], power_stdev_for_interval))
    
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

# In[]:
def compute_power_sprectra_combined (combined_power_data, sampling_freq, nperseg = 96):
    power_inst_key = 'power_inst'
    power_signal = combined_power_data[power_inst_key]/1.0e6
    #(f, S) = signal.periodogram(power_signal.flatten(), sampling_freq, scaling='density')
    (f, S) = signal.welch(power_signal.flatten(), sampling_freq, nperseg = nperseg)
    #(S, f) = plt.psd(power_signal.flatten(), Fs = sampling_freq)
    
    combined_power_data['freq'] = f
    combined_power_data['psd'] = S
    
    return combined_power_data
    