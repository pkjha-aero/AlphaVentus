#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 19:06:19 2022

@author: jha3
"""

# In[]
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
from datetime import date, datetime, timedelta, time
import pickle

# In[]
def plot_contours_instantaneous(pickle_file_name, plot_loc, qoi_plot_map, qoi_range_map, xlim=None, ylim=None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
    
    # In[] Data from 
    start_time_stamp = pickled_data_read['start_time_stamp']
    DX = pickled_data_read['DX']
    DY = pickled_data_read['DX']
    n_time_stamps = pickled_data_read['n_time_stamps']
    n_zloc = pickled_data_read['n_zloc']
    n_axial_loc = pickled_data_read['n_axial_loc']
    
    z_slices_instantaneous = pickled_data_read['z_slices_instantaneous']
    
    # In[]
    cmap_name = 'rainbow'
    
    # In[]
    for qoi in qoi_plot_map:
        qoi_unit = qoi_plot_map[qoi]
    # In[] Loop over time stamps
        for time_count, time_stamp in enumerate(z_slices_instantaneous.keys()):
            current_time_stamp = time_stamp.split('_')[1]
        
            for space_count, zloc in enumerate(z_slices_instantaneous[time_stamp]) :
                z_slice_space = z_slices_instantaneous[time_stamp][zloc]
                space = z_slice_space[qoi]
                nx = int(round(np.shape(space)[0],-2))
                ny = int(round(np.shape(space)[1],-2))
                ratio = ny/nx
                z = z_slice_space['z']
                #print(ratio)
                
                plt.figure()
                if qoi in qoi_range_map.keys():
                    qoi_range = qoi_range_map[qoi]
                    cont_levels = np.linspace(qoi_range[0], qoi_range[1], 15)
                else:
                    cont_levels = 20
                plane_cols, plane_rows = np.meshgrid(range(space.shape[1]), range(space.shape[0]))
                cont = plt.contourf(plane_cols*DX,plane_rows*DY, space, levels = cont_levels, cmap=cmap_name)
                clb = plt.colorbar(cont)
                clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
                #plt.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
                #plt.set_figwidth(5*nSpace*ratio)
                plt.xlabel(f"%s [m]"%'west_east', fontsize=14)
                plt.ylabel(f"%s [m]"%'south_north', fontsize=14)
                plt.tick_params(axis='x', labelsize=14)
                plt.tick_params(axis='y', labelsize=14)
                plt.title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
                if xlim:
                    plt.xlim([xlim[0]*DX,xlim[1]*DX])
                if ylim:
                    plt.ylim([ylim[0]*DY,ylim[1]*DY])
                plt.tight_layout() 
        
                filename = "z_slice_%s_%s_%s_%s"%(qoi, current_time_stamp, 'bottom_top', str(zloc))
                filedir = os.path.join(plot_loc, 'Instantaneous', qoi, 'bottom_top_' + str(zloc))
                if xlim or ylim:
                    filedir = os.path.join(filedir, 'zoomed')
                else:
                    filedir = os.path.join(filedir, 'full')
                os.system('mkdir -p %s'%filedir)
                plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
                plt.close()
# In[]    
def plot_contours_time_avg(pickle_file_name, plot_loc, qoi_plot_map, qoi_range_map, xlim=None, ylim=None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
    
    # In[] Data from 
    start_time_stamp = pickled_data_read['start_time_stamp']
    DX = pickled_data_read['DX']
    DY = pickled_data_read['DX']
    n_time_stamps = pickled_data_read['n_time_stamps']
    n_zloc = pickled_data_read['n_zloc']
    n_axial_loc = pickled_data_read['n_axial_loc']
    
    z_slices_time_avg = pickled_data_read['z_slices_time_avg']
    
    # In[]
    cmap_name = 'rainbow'
    
    # In[]
    for qoi in qoi_plot_map:
        qoi_unit = qoi_plot_map[qoi]
        
        time_stamp = list(z_slices_time_avg.keys())[0]
        current_time_stamp = time_stamp.split('_')[1]
    
        for space_count, zloc in enumerate(z_slices_time_avg[time_stamp]) :
            z_slice_space = z_slices_time_avg[time_stamp][zloc]
            space = z_slice_space[qoi]
            nx = int(round(np.shape(space)[0],-2))
            ny = int(round(np.shape(space)[1],-2))
            ratio = ny/nx
            z = z_slice_space['z']
            #print(ratio)
            
            plt.figure()
            if qoi in qoi_range_map.keys():
                qoi_range = qoi_range_map[qoi]
                cont_levels = np.linspace(qoi_range[0], qoi_range[1], 15)
            else:
                cont_levels = 20
            plane_cols, plane_rows = np.meshgrid(range(space.shape[1]), range(space.shape[0]))
            cont = plt.contourf(plane_cols*DX,plane_rows*DY, space, levels = cont_levels, cmap=cmap_name)
            clb = plt.colorbar(cont)
            clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
            #plt.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
            #plt.set_figwidth(5*nSpace*ratio)
            plt.xlabel(f"%s [m]"%'west_east', fontsize=14)
            plt.ylabel(f"%s [m]"%'south_north', fontsize=14)
            plt.tick_params(axis='x', labelsize=14)
            plt.tick_params(axis='y', labelsize=14)
            plt.title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
            if xlim:
                plt.xlim([xlim[0]*DX,xlim[1]*DX])
            if ylim:
                plt.ylim([ylim[0]*DY,ylim[1]*DY])
            plt.tight_layout() 
    
            filename = "z_slice_%s_%s_%s_%s"%(qoi, current_time_stamp, 'bottom_top', str(zloc))
            filedir = os.path.join(plot_loc, 'TimeAvg', qoi, 'bottom_top_' + str(zloc))
            if xlim or ylim:
                filedir = os.path.join(filedir, 'zoomed')
            else:
                filedir = os.path.join(filedir, 'full')
            os.system('mkdir -p %s'%filedir)
            plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
            plt.close()
                

# In[]:
# Plot Contours at selected time and space locations
def contourPlotSpaceTime(nc_data, qoi_data, qoi, qoi_unit, plane_header, slt, plane_rows, plane_cols, pl, plot_loc, ref_time, dt):
    desired_time = slt['Time']
    print('Desired Time : {}'.format(desired_time))
    desired_space = slt[plane_header]
    print('Desired Space ({}) : {}'.format(plane_header, desired_space))
    
    DX = nc_data.DX
    DY = nc_data.DY
    
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
   
    nSpace = np.size(desired_space)    
    nTime = np.size(desired_time) 
   
    qoi_time = []
    for time_count in range(nTime):
        qoi_space = []
        for space_count in range(nSpace):
            slice_loc = '{"%s" : %d}'%(plane_header, desired_space[space_count])
            slice_loc_dict = json.loads(slice_loc)
            arr = qoi_data.isel(Time=desired_time[time_count]).isel(slice_loc_dict)
            axes_name = arr.dims
            qoi_space.append(np.array(arr))
                    
        qoi_time.append(qoi_space)
       
    #print('qoi_time ', qoi_time)
    
        
    cmap_name = 'rainbow'
     
    if nTime==1:
        fig = plt.figure(figsize = (20, 4))
        fig.subplots_adjust(hspace=0.2, wspace=0.25)
        #fig.figure(figsize = (20, 4))
        time_inst_space = qoi_time[0]
        current_time = start_time + timedelta(seconds = desired_time[0]*dt)
        current_time_stamp = current_time.isoformat('_').split('_')[1]
        
        for space_count in range(nSpace):
            space = time_inst_space[space_count]
            nx = int(round(np.shape(space)[0],-2))
            ny = int(round(np.shape(space)[1],-2))
            ratio = ny/nx
            z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[0]).isel(bottom_top_stag = desired_space[space_count])
            z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[0]).isel(bottom_top_stag = desired_space[space_count] + 1)
            z = 0.5*(z1 + z2) 
            #print(ratio)
            ax = fig.add_subplot(nTime, nSpace, space_count+1)
            fig.set_figheight(4*nTime)
            fig.set_figwidth(5*nSpace*ratio)
            cont = ax.contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
            clb = fig.colorbar(cont, ax=ax)
            clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
            ax.set_xlabel(f"%s [m]"%axes_name[1], fontsize=14)
            ax.set_ylabel(f"%s [m]"%axes_name[0], fontsize=14)
            ax.tick_params(axis='x', labelsize=14)
            ax.tick_params(axis='y', labelsize=14)
            ax.set_title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
            ax.set_xlim([250*DX,350*DX])
            ax.set_ylim([250*DY,350*DY])
        fig.tight_layout()
    else:
        fig, ax = plt.subplots(nTime, nSpace, figsize = (20,4.5), sharex='col', sharey='row')#figsize = (20,4.5)
        fig.subplots_adjust(hspace=0.5, wspace=0.1)
    
        for time_count in range(nTime):
            time_inst_space = qoi_time[time_count]
            current_time = start_time + timedelta(seconds = desired_time[time_count]*dt)
            current_time_stamp = current_time.isoformat('_').split('_')[1]

            for space_count in range(nSpace):
                space = time_inst_space[space_count]
                nx = int(round(np.shape(space)[0],-2))
                ny = int(round(np.shape(space)[1],-2))
                ratio = ny/nx
                z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[time_count]).isel(bottom_top_stag = desired_space[space_count])
                z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[time_count]).isel(bottom_top_stag = desired_space[space_count] + 1)
                z = 0.5*(z1 + z2) 
                #print(ratio)
                cont = ax[time_count,space_count].contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
                clb = fig.colorbar(cont, ax=ax[time_count,space_count])
                clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
                fig.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
                fig.set_figwidth(5*nSpace*ratio)
                ax[time_count,space_count].set_xlabel(f"%s [m]"%axes_name[1], fontsize=14)
                ax[time_count,space_count].set_ylabel(f"%s [m]"%axes_name[0], fontsize=14)
                ax[time_count,space_count].tick_params(axis='x', labelsize=14)
                ax[time_count,space_count].tick_params(axis='y', labelsize=14)
                ax[time_count,space_count].set_title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
                ax[time_count,space_count].set_xlim([250*DX,350*DX])
                ax[time_count,space_count].set_ylim([250*DY,350*DY])
        fig.tight_layout()
    
    filename = "%s_plane_%s"%(pl, qoi)
    plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight')       
            
    
# In[]
# Plot Contours at all times and selected space locations
def contourPlotInstantaneous(nc_data, qoi_data, qoi, qoi_unit, plane_header, vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt):   
    DX = nc_data.DX
    DY = nc_data.DY
    
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
   
    nSpace = np.size(vert_loc)    
    nTime = nc_data.dims['Time']
   
    qoi_time = []
    for time_count in range(int(nTime)):
        qoi_space = []
        for space_count in range(nSpace):
            slice_loc = '{"%s" : %d}'%(plane_header, vert_loc[space_count])
            slice_loc_dict = json.loads(slice_loc)
            arr = qoi_data.isel(Time=time_count).isel(slice_loc_dict)
            axes_name = arr.dims
            qoi_space.append(np.array(arr))
                    
        qoi_time.append(qoi_space)
       
    #print('qoi_time ', qoi_time)
    
    cmap_name = 'rainbow'
     
    for time_count in range(int(nTime)):
        time_inst_space = qoi_time[time_count]
        current_time = start_time + timedelta(seconds = time_count*dt)
        current_time_stamp = current_time.isoformat('_').split('_')[1]

        for space_count in range(nSpace):
            space = time_inst_space[space_count]
            nx = int(round(np.shape(space)[0],-2))
            ny = int(round(np.shape(space)[1],-2))
            ratio = ny/nx
            z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count])
            z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count] + 1)
            z = 0.5*(z1 + z2) 
            #print(ratio)
            
            plt.figure()
            cont = plt.contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
            clb = plt.colorbar(cont)
            clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
            #plt.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
            #plt.set_figwidth(5*nSpace*ratio)
            plt.xlabel(f"%s [m]"%axes_name[1], fontsize=14)
            plt.ylabel(f"%s [m]"%axes_name[0], fontsize=14)
            plt.tick_params(axis='x', labelsize=14)
            plt.tick_params(axis='y', labelsize=14)
            plt.title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
            plt.xlim([250*DX,350*DX])
            plt.ylim([250*DY,350*DY])
            plt.tight_layout() 
   
    
            filename = "%s_plane_%s_%s_%s_%s"%(pl, qoi, current_time_stamp, plane_header, str(vert_loc[space_count]))
            plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight')       
            
# In[]:
# Plot Contours of time averaged quantities
def contourPlotTimeAvg(nc_data, qoi_data, qoi, qoi_unit, plane_header, vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt):   
    DX = nc_data.DX
    DY = nc_data.DY
    
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_').split('_')[1]
    
    nSpace = np.size(vert_loc)
    
    cmap_name = 'rainbow'

    time_count = 0
    for space_count in range(nSpace):
        space = qoi_data.isel(bottom_top = vert_loc[space_count]).mean(dim = 'Time')
        axes_name = space.dims
        nx = int(round(np.shape(space)[0],-2))
        ny = int(round(np.shape(space)[1],-2))
        ratio = ny/nx
        z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count])
        z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count] + 1)
        z = 0.5*(z1 + z2) 
        #print(ratio)
        
        plt.figure()
        cont = plt.contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
        clb = plt.colorbar(cont)
        clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
        #plt.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
        #plt.set_figwidth(5*nSpace*ratio)
        plt.xlabel(f"%s [m]"%axes_name[1], fontsize=14)
        plt.ylabel(f"%s [m]"%axes_name[0], fontsize=14)
        plt.tick_params(axis='x', labelsize=14)
        plt.tick_params(axis='y', labelsize=14)
        plt.title(f"time avg start = %s, z = %.2f m" %(start_time_stamp, z), fontsize=14)
        plt.xlim([250*DX,350*DX])
        plt.ylim([250*DY,350*DY])
        plt.tight_layout() 
   

        filename = "%s_plane_%s_time_avg_start_%s_%s_%s"%(pl, qoi, start_time_stamp, plane_header, str(vert_loc[space_count]))
        plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight') 
        
 
# In[]
# Create time-averaged line plots at selected space locations
def linePlotTimeAvg(nc_data, qoi_data, qoi, qoi_unit, we_ind, sn_ind, plot_loc, ref_time, dt):
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_').split('_')[1]
    
    time_count = 0
    z_stag = nc_data['ZTS'].isel(south_north = sn_ind).isel(west_east = we_ind).isel(Time = time_count)
    z = 0.5*(z_stag[1:] + z_stag[:-1])
    
    line_var = qoi_data.isel(south_north = sn_ind).isel(west_east = we_ind).mean(dim = 'Time')
    
    plt.figure()
    plt.plot(line_var, z)
    plt.xlabel(f"%s [%s]"%(qoi, qoi_unit), fontsize=14)
    plt.ylabel(f"Height [m]"%z, fontsize=14)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title(f"time avg start = %s, we = %d, sn = %d m" %(start_time_stamp, we_ind, sn_ind), fontsize=14)
    #plt.xlim([250*DX,350*DX])
    plt.ylim([0, 150])
    plt.tight_layout() 


    filename = "vertline_%s_we_%s_sn_%s_time_avg_start_%s"%(qoi, str(we_ind), str(sn_ind),  start_time_stamp)
    plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight')
    
    