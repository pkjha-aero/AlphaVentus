#!/usr/bin/env python
# coding: utf-8

# In[]:
# Import Modules
from data_processing import *
from combine_data import *
from plotting import *

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
                          'part15': '2010-05-16_01:40:10'
    }
#'''
    
# In[]:
z_plane_locs = [20]

vert_line_locs = [1] #D

dt = 10 # sec

frac_time = 0.02 # Fraction of time series to use
frac_time = 1.00 # Fraction of time series to use

# In[]:
case_name = 'MesoMicro1_CPM'
processed_results_loc_base = '/Users/jha3/Downloads/AlphaVentusProcessed_PowerOnly'
processed_results_loc = os.path.join(processed_results_loc_base, case_name)
proc_data_loc_for_combined = os.path.join(processed_results_loc, 'combined_procfiles')
 
# pickle file for combined
pickle_file_name_for_combined = '{}_pickled.pkl'.format('combined_data')
pickle_file_loc_for_combined = os.path.join(proc_data_loc_for_combined, pickle_file_name_for_combined)
pickled_data_combined = {}

# In[]
for interval, tsout_file_stamp in interval_tsoutfile_map.items():   
    # Time stamps for tsoutfile and interval
    tsout_time = tsout_file_stamp.split('_')
    start_time = datetime.fromisoformat(tsout_time[-2]+ '_' + tsout_time[-1]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_') #.split('_')[1]
    
    # Processed data location for this interval
    proc_data_loc_for_interval = os.path.join(processed_results_loc, '{}_procfiles'.format(interval))
   
    # pickle file for interval
    pickle_file_name_for_interval = '{}_pickled_{}.pkl'.format(interval, start_time_stamp)
    pickle_file_loc_for_interval = os.path.join(proc_data_loc_for_interval, pickle_file_name_for_interval)
    
    # Open pickled data for the interval
    with open(pickle_file_loc_for_interval, 'rb') as pickle_file_handle:
        pickled_data_for_interval = pickle.load(pickle_file_handle)
    
    # Combine the picked data for intervals
    pickled_data_combined = combine_power_for_intervals(pickled_data_combined, pickled_data_for_interval)
    

# In[]
# Write the combined pickled data
os.system('mkdir -p %s'%proc_data_loc_for_combined)
with open(pickle_file_loc_for_combined, 'wb') as pickle_file_handle:
    pickle.dump(pickled_data_combined, pickle_file_handle)

# In[]    
plot_power_inst (pickle_file_loc_for_combined, proc_data_loc_for_combined, case_name, dt)
plot_power_avg (pickle_file_loc_for_combined, proc_data_loc_for_combined, case_name)
plot_power_histogram (pickle_file_loc_for_combined, proc_data_loc_for_combined, case_name, num_bins = 30)
    
# In[]
# Open the combined pickled data saved above
'''
with open(pickle_file_loc_for_combined, 'rb') as pickle_file_handle:
    pickled_data_combined_read = pickle.load(pickle_file_handle)
'''
# In[]
sim_end_time = timer()
print('Total computing time: {} s'.format(sim_end_time - sim_start_time))

# In[]

# Loop over QoI of interest to create plots
for (qoi, qoi_unit) in zip(qoi_list, qoi_units):
    #qoi_data, qoi_dim_names, qoi_dim = getQoI(tsout_domain6, qoi)
    # In[]    
    '''
    # Plot contours at the desired times and spatial locations
    slt = sampling_loc_time(qoi_dim_names, qoi_dim)
    for plane_ind in range(0,1):
        pl = plane[plane_ind]
        if pl=='xy':
            plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
            contourPlotSpaceTime(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1],slt,plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
        elif pl=='yz':
            plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[2]], wrf_domain6[qoi_dim_names[1]])
            contourPlotSpaceTime(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[3],slt,plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
        elif pl=='xz':
            plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[1]])
            contourPlotSpaceTime(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[2],slt,plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    
    # In[]  
    '''
    pl = 'xy'
    vert_loc = [20, 25]
    plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
    contourPlotInstantaneous(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1], vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    
    # In[]  
    '''
    pl = 'xy'
    vert_loc = [20]
    plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
    contourPlotTimeAvg(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1], vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    
    # In[]
    '''
    [we_ind, sn_ind] = [300, 300]
    linePlotTimeAvg(wrf_domain6, qoi_data, qoi, qoi_unit, we_ind, sn_ind, plot_loc, ref_time, dt)
    '''
    
    # In[]
    
# In[]
dummy = 0
