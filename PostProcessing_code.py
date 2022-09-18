#!/usr/bin/env python
# coding: utf-8

# In[]:
# Import Modules
#from data_processing import *
#from plotting import *

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

plot_loc = '/Users/jha3/Downloads/plots'

# In[]:
# reference time (based on NetCDF file name) and sampling time interval
ref_time = netCDF_file_loc.split('_')
dt = 10 # sec
pickle_file_name = 'pickle_' + ref_time[2] + '_' + ref_time[3]

# In[]:
# Variables of interest
plane = ['xy','yz','xz']

qoi_list = ['UTS', 'TKETS', 'T12TS']
#qoi_list = ['UTS']
qoi_units = ['m/s', 'm2/s2', 'm2/s2']

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
    [we_ind, sn_ind] = [300, 300]
    linePlotTimeAvg(wrf_domain6, qoi_data, qoi, qoi_unit, we_ind, sn_ind, plot_loc, ref_time, dt)
    
    # In[]  
    dummy = 0