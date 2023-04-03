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
    
    