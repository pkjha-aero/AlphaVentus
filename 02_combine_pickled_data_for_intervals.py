#!/usr/bin/env python
# coding: utf-8

# In[]:
# Import Modules
from data_processing_base import *
from data_processing_power import *
from data_processing_slice import *
from plotting_power import *
from plotting_slice import *
from plotting_line  import *

# In[]
import os
import sys
#import wrf
import numpy as np
import math
import xarray as xr
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
import pickle
from datetime import date, datetime, timedelta, time

from timeit import default_timer as timer

plt.style.use('seaborn-white')

# In[]
sim_start_time = timer()

# In[]:
# Variables of interest

qoi_from_tsout_file = ['UTS','VTS', 'WTS', 'TKETS']

qoi_plot_map = {'UMAG': 'm/s',
                'TKE_RES': 'm2/s2',
                'TKE_SGS': 'm2/s2',
                'TKE_TOT': 'm2/s2',
                'T23': 'm2/s2',
                'T13': 'm2/s2',
                'W'   : 'm/s'
               }
qoi_plot_avg_map = {'UMAG_AVG': 'm/s',
                    'TKE_RES_AVG': 'm2/s2',
                    'TKE_SGS_AVG': 'm2/s2',
                    'TKE_TOT_AVG': 'm2/s2',
                    'W_AVG'   : 'm/s'
                   }

qoi_range_map = {'UMAG'    : [8.0, 20.0],
                 'UMAG_AVG': [8.0, 20.0],
                 'TKE_RES' : [0.0, 10.0],
                 'TKE_RES_AVG' : [0.0, 10.0],
                 'TKE_SGS' : [0.0, 0.6],
                 'TKE_SGS_AVG' : [0.0, 0.6],
                 'W'       : [-1.35, 1.35],
                 'W_AVG'   : [-1.35, 1.35]
                }

# In[]:
interval_tsoutfile_map = {'part05': '2010-05-16_00:00:10'
        }
#'''    
interval_tsoutfile_map = {'part05': '2010-05-16_00:00:10',
                          'part06': '2010-05-16_00:10:10',
                          'part07': '2010-05-16_00:20:10',
                          'part08': '2010-05-16_00:30:10',
                          'part09': '2010-05-16_00:40:10',
                          'part10': '2010-05-16_00:50:10',
                          'part11': '2010-05-16_01:00:10',
                          'part12': '2010-05-16_01:10:10',
                          'part13': '2010-05-16_01:20:10',
                          'part14': '2010-05-16_01:30:10',
                          'part15': '2010-05-16_01:40:10',
                          'part16': '2010-05-16_01:50:10'
    }
#'''
    
# In[]:
z_plane_locs = [20]

vert_line_locs = [1] #D

dt = 10 # sec

#frac_time = 0.02 # Fraction of time series to use
#frac_time = 1.00 # Fraction of time series to use

# In[]:
# Flags etc.
combine_slice_for_intervals = True
plot_slice_combined_intervals = True

# In[]
combine_power_for_intervals = True
scale_power_for_intervals = True
plot_power_combined_intervals = True

compute_power_stdev = True
plot_power_stdev = True

compute_power_pdf = True

compute_power_spectra = True
plot_power_spectra = True


# In[]:
case_name = 'MesoMicro1_CPM'

# Path where processed data for intervals and combined data are located
processed_slice_base_loc = '/Users/jha3/Downloads/AlphaVentus_Slice'
processed_slice_case_loc = os.path.join(processed_slice_base_loc, case_name)

processed_power_base_loc = '/Users/jha3/Downloads/AlphaVentus_Power'
processed_power_case_loc = os.path.join(processed_power_base_loc, case_name)

# In[]:
if combine_slice_for_intervals or plot_slice_combined_intervals:
    combined_slice_data = {}
    combined_slice_data_loc = os.path.join(processed_slice_case_loc, 'combined_procfiles')
    combined_slice_pickle_file = os.path.join(combined_slice_data_loc, '{}_combined_slice.pkl'.format(case_name))

if combine_power_for_intervals or plot_power_combined_intervals:
    combined_power_data = {}
    combined_power_data_loc = os.path.join(processed_power_case_loc, 'combined_procfiles')
    combined_power_pickle_file = os.path.join(combined_power_data_loc, '{}_combined_power.pkl'.format(case_name))

# In[]
for interval, tsout_file_stamp in interval_tsoutfile_map.items():   
    # Time stamps for tsoutfile and interval
    tsout_time = tsout_file_stamp.split('_')
    start_time = datetime.fromisoformat(tsout_time[-2]+ '_' + tsout_time[-1]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_') #.split('_')[1]
   
    # In[]:
    if combine_slice_for_intervals:
        # Processed slice data location for this interval
        processed_slice_interval_loc = os.path.join(processed_slice_case_loc, '{}_procfiles'.format(interval))
        
        # Pickle file name for processed slice data
        pickled_slice_interval_file_name = '{}_slice_{}.pkl'.format(interval, start_time_stamp)
        pickled_slice_interval_file = os.path.join(processed_slice_interval_loc, pickled_slice_interval_file_name)
    
    # In[]:
    if combine_power_for_intervals:
        # Processed slice data location for this interval
        processed_power_interval_loc = os.path.join(processed_power_case_loc, '{}_procfiles'.format(interval))
        
        # Pickle file name for processed power data
        pickled_power_interval_file_name = '{}_power_{}.pkl'.format(interval, start_time_stamp)
        pickled_power_interval_file = os.path.join(processed_power_interval_loc, pickled_power_interval_file_name)
    
        # Open pickled data for the interval
        with open(pickled_power_interval_file, 'rb') as pickled_power_interval_file_handle:
            pickled_power_data_interval = pickle.load(pickled_power_interval_file_handle)
            
        # Scale power data if desired
        if scale_power_for_intervals:
            pickled_power_data_interval = scale_power(pickled_power_data_interval, 4.0)
    
        # Combine the picked data for intervals
        combined_power_data = combine_power_for_intervals_of_case(combined_power_data, pickled_power_data_interval, compute_power_stdev)
 
# In[]
if compute_power_pdf:
    combined_power_data = compute_power_pdf_combined (combined_power_data, num_bins=30)
    
# In[]:
if compute_power_spectra:
    combined_power_data = compute_power_sprectra_combined(combined_power_data, sampling_freq = 1.0/dt, nperseg = 96)
    
# In[]
# Write the combined pickled data
if combine_power_for_intervals:
    os.system('mkdir -p %s'%combined_power_data_loc)
    with open(combined_power_pickle_file, 'wb') as combined_power_pickle_file_handle:
        pickle.dump(combined_power_data, combined_power_pickle_file_handle)
# In[]
if plot_power_combined_intervals:
    '''
    plot_power_inst(combined_power_pickle_file, combined_power_data_loc, case_name, dt)
    plot_power_inst(combined_power_pickle_file, combined_power_data_loc, case_name, dt, [2, 8])
    plot_power_pdf(combined_power_pickle_file, combined_power_data_loc, case_name, 30)
    plot_power_pdf(combined_power_pickle_file, combined_power_data_loc, case_name, 30, [2, 8], [0, 0.5])
    '''
    plot_power_inst_pdf (combined_power_pickle_file, combined_power_data_loc, case_name, dt, 30)
    plot_power_inst_pdf (combined_power_pickle_file, combined_power_data_loc, case_name, dt, 30, [2.5, 5.5], [0, 20.0])
    plot_power_inst_pdf (combined_power_pickle_file, combined_power_data_loc, case_name, dt, 30, [2.5, 5.5], None)
    plot_power_inst_pdf (combined_power_pickle_file, combined_power_data_loc, case_name, dt, 30, None, [0, 20.0])
    
    plot_power_avg(combined_power_pickle_file, combined_power_data_loc, case_name)
    plot_power_avg(combined_power_pickle_file, combined_power_data_loc, case_name, [3.5, 5.0])
    
    if compute_power_stdev and plot_power_stdev:
        plot_power_stdev_combined(combined_power_pickle_file, combined_power_data_loc, case_name)
        
    if compute_power_spectra and plot_power_spectra:
        plot_power_spectra_combined(combined_power_pickle_file, combined_power_data_loc, case_name)
    
# In[]
sim_end_time = timer()
print('Total computing time: {} s'.format(sim_end_time - sim_start_time))

# In[]
dummy = 0
