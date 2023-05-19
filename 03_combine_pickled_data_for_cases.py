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
combine_slice_for_cases = False
combine_power_for_cases = True

plot_power_for_cases = True

# In[]:
case_names = ['MesoMicro1_CPM', 'MesoMicro1_NO_CPM', 'MesoMicro2_CPM', 'MesoMicro2_NO_CPM']
case_colors = {'MesoMicro1_CPM': 'r', 
               'MesoMicro1_NO_CPM': 'b',
               'MesoMicro2_CPM': 'g',
               'MesoMicro2_NO_CPM': 'm'}

# In[]
# Path where processed data for intervals and combined data are located
processed_power_base_loc = '/Users/jha3/Downloads/AlphaVentus_Power'


# In[]:
if combine_power_for_cases:
    combined_power_data_all_cases = {}
    combined_power_data_loc = os.path.join(processed_power_base_loc, 'combined_procfiles_all_cases')
    combined_power_pickle_file = os.path.join(combined_power_data_loc, 'combined_power_all_cases.pkl')

# In[]:
for case_name in case_names:
    if combine_power_for_cases:
        processed_power_case_loc = os.path.join(processed_power_base_loc, case_name)
        combined_power_case_pickle_loc = os.path.join(processed_power_case_loc, 'combined_procfiles')
        combined_power_case_pickle_file = os.path.join(combined_power_case_pickle_loc, '{}_combined_power.pkl'.format(case_name))
        
        with open(combined_power_case_pickle_file, 'rb') as combined_power_case_pickle_file_handle:
            combined_power_data_case = pickle.load(combined_power_case_pickle_file_handle)
        combined_power_data_all_cases[case_name] = combined_power_data_case
    
# In[]
# Write the combined pickled data
if combine_power_for_cases:
    os.system('mkdir -p %s'%combined_power_data_loc)
    with open(combined_power_pickle_file, 'wb') as combined_power_pickle_file_handle:
        pickle.dump(combined_power_data_all_cases, combined_power_pickle_file_handle)
# In[]
if plot_power_for_cases:
    plot_power_avg_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, 'All')
    plot_power_avg_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, 'All', [3.9, 4.6])
    
    plot_power_stdev_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, 'All')
    
    plot_power_pdf_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, case_names, 'All')
    plot_power_pdf_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, ['MesoMicro1_CPM', 'MesoMicro2_CPM'], 'CPM')
    plot_power_pdf_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, ['MesoMicro1_NO_CPM', 'MesoMicro2_NO_CPM'], 'No_CPM')
    
    plot_power_spectra_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, case_names, 'All')
    plot_power_spectra_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, ['MesoMicro1_CPM', 'MesoMicro2_CPM'], 'CPM')
    plot_power_spectra_all_cases(combined_power_pickle_file, combined_power_data_loc, case_colors, ['MesoMicro1_NO_CPM', 'MesoMicro2_NO_CPM'], 'No_CPM')

# In[]
sim_end_time = timer()
print('Total computing time: {} s'.format(sim_end_time - sim_start_time))

# In[]
dummy = 0
