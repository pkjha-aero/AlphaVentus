#!/usr/bin/env python
# coding: utf-8

# In[]:
# Import Modules
from data_processing import *
from plotting import *

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
plt.style.use('seaborn-white')

# In[]
# Load NetCDF file containing tsout data and define plot location
#scratch_folder_loc = '/p/lustre2/nainap/'
#netCDF_file_loc = os.path.join(scratch_folder_loc,'From_Pankaj', 'NetCDF_files','wrfout_GAD_FINO_SecondResult' )
netCDF_file_loc = '/Users/jha3/Downloads/tsout_d06_2010-05-16_00:00:10'
wrf_domain6 = xr.open_dataset(netCDF_file_loc)

proc_data_loc = '/Users/jha3/Downloads/proc_data'

# In[]:
# reference time (based on NetCDF file name) and sampling time interval
ref_time = netCDF_file_loc.split('_')
dt = 10 # sec

# In[]:
# Variables of interest
plane = ['xy','yz','xz']

qoi_units_map = {'UTS': 'm/s',
                 'VTS': 'm/s',
                 'WTS': 'm/s'
                }
qoi_list = list(qoi_units_map.keys())
#qoi_list = ['UTS']
qoi_units = list(qoi_units_map.values())

qoi_plot_map = {'UMAG': 'm/s',
                #'TKE' : 'm2/s2',
                'W'   : 'm/s'
               }

z_plane_locs = [20, 25]

vert_line_locs = [1] #D

frac_time = 0.05 # Fraction of tiem series to use

'''
# In[]
pickled_data, pickle_file = create_instantaneous_data (wrf_domain6, qoi_units_map, z_plane_locs, vert_line_locs, ref_time, dt, frac_time)
#create_instantaneous_data (wrf_domain6, qoi_units_map, z_plane_locs, vert_line_locs, ref_time, dt, frac_time)

with open(os.path.join(proc_data_loc, pickle_file), 'wb') as pickle_file_handle:
    pickle.dump(pickled_data, pickle_file_handle)
'''    
# In[]
pickle_file = 'pickled_2010-05-16_00:00:00.pkl'
plot_contours_instantaneous(os.path.join(proc_data_loc, pickle_file), proc_data_loc, qoi_plot_map)

# In[]
# Loop over QoI of interest to create plots
for (qoi, qoi_unit) in zip(qoi_list, qoi_units):
    qoi_data, qoi_dim_names, qoi_dim = getQoI(wrf_domain6, qoi)

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