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
def get_z_slices_instantaneous(nc_data, qoi_from_tsout_file, z_plane_locs, ref_time, dt, start_time, frac_time):
    z_slices_instantaneous = {}
    
    n_time_stamps = nc_data.dims['Time']
    n_zloc        = np.size(z_plane_locs)
    
    z_slice_data = {}
    
    # In[]
    # In[]
    qoi_list = ['UTS', 'VTS', 'WTS', 'T13TS', 'T23TS']
    #qoi_list = ['UTS', 'VTS']
    #qoi_from_tsout_file = qoi_list
    z_slices_time_avg = get_z_slices_time_averaged(nc_data, qoi_list, z_plane_locs, ref_time, dt, start_time, frac_time)
   
    
    # In[] Read all the relevant data from the NetCDF file
    for qoi in qoi_from_tsout_file:
        qoi_data = nc_data[qoi]
        if qoi in ['WTS', 'ZTS', 'T13TS', 'T23TS']:
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
            z_slice_space['U'] = z_slice_data['UTS'][time_count][space_count]
            z_slice_space['V'] = z_slice_data['VTS'][time_count][space_count]
            z_slice_space['W'] = z_slice_data['WTS'][time_count][space_count]
            z_slice_space['UMAG'] = z_slice_data['UMAG'][time_count][space_count]
            z_slice_space['TKE_SGS'] = z_slice_data['TKETS'][time_count][space_count]
            z_slice_space['T23'] = z_slice_data['T23TS'][time_count][space_count]
            z_slice_space['T13'] = z_slice_data['T13TS'][time_count][space_count]
            
            z_slice_space['UP'] = z_slice_space['U'] - z_slices_time_avg[start_time.isoformat('_')][z_plane_locs[space_count]]['U_AVG']
            z_slice_space['VP'] = z_slice_space['V'] - z_slices_time_avg[start_time.isoformat('_')][z_plane_locs[space_count]]['V_AVG']
            z_slice_space['WP'] = z_slice_space['W'] - z_slices_time_avg[start_time.isoformat('_')][z_plane_locs[space_count]]['W_AVG']
            z_slice_space['TKE_RES'] = 0.5*(np.square(z_slice_space['UP']) + np.square(z_slice_space['VP']) + np.square(z_slice_space['WP']))
            z_slice_space['TKE_TOT'] = z_slice_space['TKE_RES'] + z_slice_space['TKE_SGS']
            
            z_slice_time[z_plane_locs[space_count]] = z_slice_space
        
        z_slices_instantaneous[current_time_stamp] = z_slice_time
            
    # In[]
    for qoi in ['UMAG', 'TKE_RES', 'TKE_SGS', 'TKE_TOT']:
    #for qoi in ['UMAG']:
        timestamp_key = list(z_slices_instantaneous)[0]
        zloc_key = list(z_slices_instantaneous[timestamp_key])[0]
        space_dim = z_slices_instantaneous[timestamp_key][zloc_key]['UMAG'].shape
        space_avg_zloc = np.zeros((n_zloc, space_dim[0], space_dim[1]))
        
        for time_count, time_stamp in enumerate(z_slices_instantaneous.keys()):
            current_time_stamp = time_stamp.split('_')[1]
        
            for space_count, zloc in enumerate(z_slices_instantaneous[time_stamp]): 
                z_slice_space = z_slices_instantaneous[time_stamp][zloc]
                space = z_slice_space[qoi]
                if time_count == 0:
                    space_avg_zloc[space_count] = space
                else:
                    space_avg_zloc[space_count] += space
         
        for space_count in range(space_avg_zloc.shape[0]):
            space_avg = space_avg_zloc[space_count]/len(z_slices_instantaneous.keys())
            z_slices_time_avg[start_time.isoformat('_')][z_plane_locs[space_count]][qoi + '_AVG'] = space_avg
            
    # In[]
    return z_slices_instantaneous, z_slices_time_avg

# In[]
def get_z_slices_time_averaged(nc_data, qoi_list, z_plane_locs, ref_time, dt, start_time, frac_time):
    z_slices_time_averaged = {}
    
    n_time_stamps = nc_data.dims['Time']
    n_zloc        = np.size(z_plane_locs)
    
    z_slice_data = {}
    time_count = 0
    
    # In[] Read all the relevant data from the NetCDF file
    for qoi in qoi_list:
        qoi_data = nc_data[qoi]
        if qoi in ['WTS', 'ZTS', 'T13TS', 'T23TS']:
            bottom_top_dim = 'bottom_top_stag'
        else:
            bottom_top_dim = 'bottom_top'
        
        z_slice_data[qoi] = []
        z_slice_avg = []
        for space_count in range(n_zloc):
            slice_loc = '{"%s" : %d}'%(bottom_top_dim, z_plane_locs[space_count])
            slice_loc_dict = json.loads(slice_loc)
            arr = qoi_data.isel(slice_loc_dict).mean(dim = 'Time')
            z_slice_avg.append(np.array(arr))
        z_slice_data[qoi].append(z_slice_avg)
        
    # In[] Rearrange data in dict format
    current_time_stamp = start_time.isoformat('_')
    z_slice_time = {}
    
    for space_count in range(n_zloc):
        z_slice_space = {}
        z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = z_plane_locs[space_count])
        z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = z_plane_locs[space_count] + 1)
        z_slice_space['z'] = 0.5*np.float((z1 + z2))
        z_slice_space['U_AVG'] = z_slice_data['UTS'][time_count][space_count]
        z_slice_space['V_AVG'] = z_slice_data['VTS'][time_count][space_count]
        z_slice_space['W_AVG'] = z_slice_data['WTS'][time_count][space_count]
        z_slice_space['T23_AVG'] = z_slice_data['T23TS'][time_count][space_count]
        z_slice_space['T13_AVG'] = z_slice_data['T13TS'][time_count][space_count]

        z_slice_time[z_plane_locs[space_count]] = z_slice_space
    
    z_slices_time_averaged[current_time_stamp] = z_slice_time
            
    # In[]    
    return z_slices_time_averaged

# In[]
def create_slice_data (nc_data, qoi_from_tsout_file, z_plane_locs, vert_line_locs, ref_time, dt, start_time, frac_time):
    pickled_data = create_base_data(nc_data, start_time)
    
    n_zloc        = np.size(z_plane_locs)
    pickled_data.update({'n_zloc': n_zloc})
    
    pickled_data['z_slices_instantaneous'], pickled_data['z_slices_time_avg'] = get_z_slices_instantaneous(nc_data, qoi_from_tsout_file, z_plane_locs, ref_time, dt, start_time, frac_time)
    
    return pickled_data
   