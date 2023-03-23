#!/usr/bin/env python
# coding: utf-8

# In[]:
# Import Modules
from data_processing import *
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

qoi_plot_map = {'UMAG': 'm/s'}

qoi_plot_avg_map = {'UMAG_AVG': 'm/s',
                    'TKE_RES_AVG': 'm2/s2',
                    'TKE_SGS_AVG': 'm2/s2',
                    'TKE_TOT_AVG': 'm2/s2',
                    'W_AVG'   : 'm/s'
                   }

qoi_plot_avg_map = {'UMAG_AVG': 'm/s'}

qoi_range_map = {'UMAG'    : [8.0, 15.0],
                 'UMAG_AVG': [8.0, 15.0],
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
'''    
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
'''
    
# In[]:
z_plane_locs = [20]

vert_line_locs = [1] #D

dt = 10 # sec

frac_time = 0.05 # Fraction of time series to use
#frac_time = 1.00 # Fraction of time series to use

# In[]:
case_name = 'MesoMicro1_CPM'

WRF_result_loc_base ='/Users/jha3/Downloads/AlphaVentusSimOutput'
WRF_result_files_loc = os.path.join(WRF_result_loc_base, case_name)

processed_results_loc_base = '/Users/jha3/Downloads/AlphaVentusProcessed_Test'
processed_results_loc = os.path.join(processed_results_loc_base, case_name)

# In[]
for interval, tsout_file_stamp in interval_tsoutfile_map.items():
    interval_dir = '{}_outfiles'.format(interval)
    tsoutfile_name = 'tsout_d06_{}'.format(tsout_file_stamp)
    netCDF_file_loc = os.path.join(WRF_result_files_loc, interval_dir, tsoutfile_name)
    
    # Time stamps for tsoutfile and interval
    tsout_time = netCDF_file_loc.split('_')
    start_time = datetime.fromisoformat(tsout_time[-2]+ '_' + tsout_time[-1]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_') #.split('_')[1]
    
    # Processed data location for this interval
    proc_data_loc = os.path.join(processed_results_loc, '{}_procfiles'.format(interval))
    
    # pickle file name
    pickle_file_name = '{}_pickled_{}.pkl'.format(interval, start_time_stamp)
    pickle_file_loc = os.path.join(proc_data_loc, pickle_file_name)
    #'''
    # Open tsout NetCDF data
    tsout_domain6 = xr.open_dataset(netCDF_file_loc)
    
    # Create pickled data
    pickled_data = create_data (tsout_domain6, qoi_from_tsout_file, z_plane_locs, vert_line_locs, tsout_time, dt, start_time, frac_time)
    
    # Write picked data
    os.system('mkdir -p %s'%proc_data_loc)
    with open(pickle_file_loc, 'wb') as pickle_file_handle:
        pickle.dump(pickled_data, pickle_file_handle)
    #'''
    '''
    # Plot Line Plot of instantaneous and avg power
    plot_power_inst (pickle_file_loc, proc_data_loc, case_name, dt)
    plot_power_avg (pickle_file_loc, proc_data_loc, case_name)
    '''
    # Plot contour plots of instantaneous data
    plot_contours_instantaneous(pickle_file_loc, proc_data_loc, qoi_plot_map, qoi_range_map, [250, 350], [250, 350])
    plot_contours_instantaneous(pickle_file_loc, proc_data_loc, qoi_plot_map, qoi_range_map)
    
    # Plot contour plots of averaged data
    plot_contours_time_avg(pickle_file_loc, proc_data_loc, qoi_plot_avg_map, qoi_range_map, [250, 350], [250, 350])
    plot_contours_time_avg(pickle_file_loc, proc_data_loc, qoi_plot_avg_map, qoi_range_map)

# In[]
sim_end_time = timer()
print('Total computing time: {} s'.format(sim_end_time - sim_start_time))

    
# In[]
dummy = 0
