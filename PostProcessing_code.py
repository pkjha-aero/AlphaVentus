#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:
from helper_functions import *

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
plt.style.use('seaborn-white')


# # Load NetCDF From Cluster

# In[2]:


#scratch_folder_loc = '/p/lustre2/nainap/'
#netCDF_file_loc = os.path.join(scratch_folder_loc,'From_Pankaj', 'NetCDF_files','wrfout_GAD_FINO_SecondResult' )


# In[3]:


netCDF_file_loc = '/Users/jha3/Downloads/tsout_d06_2010-05-16_00:00:10'


# In[4]:


wrf_domain6 = xr.open_dataset(netCDF_file_loc)


# In[5]:


#wrf_domain6


# In[6]:


#wrf_domain6['UTS']


# # Set-Up Variables for Contour Plots


# In[]:
ref_time = netCDF_file_loc.split('_')
dt = 10 # sec


# In[10]:


plane = ['xy','yz','xz']
qoi_list = ['UTS', 'TKETS', 'T12TS']
#qoi_list = ['UTS']
qoi_units = ['m/s', 'm2/s2', 'm2/s2']

plot_loc = '/Users/jha3/Downloads/plots'

for (qoi, qoi_unit) in zip(qoi_list, qoi_units):
    qoi_data, qoi_dim_names, qoi_dim = getQoI(wrf_domain6, qoi)

        
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
    '''
    pl = 'xy'
    vert_loc = [20, 25]
    plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
    contourPlotInstantaneous(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1], vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    '''
    pl = 'xy'
    vert_loc = [20]
    plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
    contourPlotTimeAvg(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1], vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    [we_ind, sn_ind] = [300, 300]
    linePlotTimeAvg(wrf_domain6, qoi_data, qoi, qoi_unit, we_ind, sn_ind, plot_loc, ref_time, dt)
    
    dummy = 0