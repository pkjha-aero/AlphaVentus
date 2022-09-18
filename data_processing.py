#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 03:13:43 2022

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
def getQoI (netCDF_data, qoi):
    qoi_data = netCDF_data[qoi]
    qoi_dim_names = qoi_data.dims
    qoi_dim = qoi_data.shape
    
    #shape format in : [time, z, y, x]
    print('QoI: {}, \nDimension names: {}, \nDimensions: {} '.format(qoi, qoi_dim_names, qoi_dim))
    return qoi_data, qoi_dim_names, qoi_dim

# In[]:
def sampling_loc_time(qoi_dim_names, qoi_dim):
    [time_dim, bottom_top_dim, south_north_dim, west_east_dim] = qoi_dim
    
    desired_time = [0, time_dim-1]
    #desired_time = list(np.arange(0, time_dim, 1))
    #hub_height = int(wrf_domain6['HUB_HEIGHT'].isel(Time=0)[0])
    hub_height = 100
    #bottom_top = [0.05*hub_height, 0.1*hub_height, 0.5*hub_height,hub_height] #bottom_top_locs
    bottom_top = [15,20,25]
    south_north = [int(west_east_dim*0.5), int(west_east_dim*0.75)] #south_north_locs
    west_east = [int(bottom_top_dim*0.25), int(south_north_dim*0.5), int(south_north_dim*0.75)] #west_east_locs

    #sampling location dictionary
    sampling_loc_time = '{"%s" : %s, "%s" : %s, "%s" : %s, "%s" : %s}'%(qoi_dim_names[0],desired_time, qoi_dim_names[1],bottom_top, qoi_dim_names[2], south_north, qoi_dim_names[3],west_east)
    slt = json.loads(sampling_loc_time)
    
    return slt

# In[]
def get_z_slices_instantaneous(nc_data, qoi_units_map, z_plane_locs, ref_time, dt, start_time, frac_time):
    z_slices_instantaneous = {}
    
    n_time_stamps = nc_data.dims['Time']
    n_zloc        = np.size(z_plane_locs)
    
    z_slice_data = {}
    
    # In[] Read all the relevant data from the 
    for qoi in qoi_units_map.keys():
        qoi_data = nc_data[qoi]
        if qoi == 'WTS' or qoi == 'ZTS':
            bottom_top_dim = 'bottom_top_stag'
        else:
            bottom_top_dim = 'bottom_top'
        
        z_slice_data[qoi] = []
        for time_count in range(int(frac_time*n_time_stamps)):
            z_slice_t = []
            for space_count in range(n_zloc):
                slice_loc = '{"%s" : %d}'%(bottom_top_dim, z_plane_locs[space_count])
                slice_loc_dict = json.loads(slice_loc)
                arr = qoi_data.isel(Time=time_count).isel(slice_loc_dict)
                axes_name = arr.dims
                z_slice_t.append(np.array(arr))
            z_slice_data[qoi].append(z_slice_t)
     
    # In[] Get U magnitude
    z_slice_data['UMAG'] = []        
    for time_count in range(int(frac_time*n_time_stamps)):
        z_slice_t = []
        for space_count in range(n_zloc):
            z_slice_u = z_slice_data['UTS'][time_count][space_count]
            z_slice_v = z_slice_data['VTS'][time_count][space_count]
            z_slice_umag = np.sqrt(np.square(z_slice_u) + np.square(z_slice_v))
            z_slice_t.append(z_slice_umag)
        z_slice_data['UMAG'].append(z_slice_t)
        
    # In[] Rearrange data in dict format
    for time_count in range(int(frac_time*n_time_stamps)):
        current_time = start_time + timedelta(seconds = time_count*dt)
        current_time_stamp = current_time.isoformat('_')
        z_slice_time = {}
        
        for space_count in range(n_zloc):
            z_slice_space = {}
            z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = z_plane_locs[space_count])
            z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = z_plane_locs[space_count] + 1)
            z_slice_space['z'] = 0.5*np.float((z1 + z2))
            z_slice_space['W'] = z_slice_data['WTS'][time_count][space_count]
            z_slice_space['UMAG'] = z_slice_data['UMAG'][time_count][space_count]

            z_slice_time[z_plane_locs[space_count]] = z_slice_space
        
        z_slices_instantaneous[current_time_stamp] = z_slice_time
            
    # In[]    
    return z_slices_instantaneous

# In[]
def create_instantaneous_data (nc_data, qoi_units_map, z_plane_locs, vert_line_locs, ref_time, dt, frac_time):
    pickled_data = {}
    
    # Start time stamp to identify the beginning of a data set 
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_') #.split('_')[1]
    pickled_data['start_time_stamp'] = start_time_stamp
    pickle_filename = f"pickled_%s.pkl"%start_time_stamp
    
    # Extract grid resolution
    DX = nc_data.DX
    DY = nc_data.DY
    pickled_data.update({'DX': DX, 'DY': DY})
    
    # Find number of time stamps, z-locations, and axial locations for sampling
    n_time_stamps = nc_data.dims['Time']
    n_zloc        = np.size(z_plane_locs)
    n_axial_loc   = np.size(vert_line_locs)
    pickled_data.update({'n_time_stamps': n_time_stamps, 'n_zloc': n_zloc, 'n_axial_loc': n_axial_loc})
    
    pickled_data['z_slices_instantaneous'] = get_z_slices_instantaneous(nc_data, qoi_units_map, z_plane_locs, ref_time, dt, start_time, frac_time)
    
    return pickled_data, pickle_filename
    
    