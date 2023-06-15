# Import Modules
# In[]:

import sys
import os
import os.path as path
import xarray as xr
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# In[]:
sys.path.insert(0, os.path.dirname(os.getcwd()))
from data_processing_slice import *
from plotting_slice import *

# In[]:
def compute_indices_weights_for_slices (case_dir_map, hub_height):
    for case in case_dir_map.keys():
        dz = case_dir_map[case]['dz']

        z_ind1 = hub_height // int(dz)
        z_ind2 = z_ind1 + 1
        z1 = int(dz*z_ind1)
        z2 = int(dz*z_ind2)
        w1 = 1.0 - (hub_height - z1)/dz
        w2 = 1.0 - (z2 - hub_height)/dz
        print('z_indices: [{}, {}], z: [{}, {}], weights: [{}, {}]'.format(
            z_ind1, z_ind2, z1, z2, w1, w2))

        case_dir_map[case]['z_indices'] = [z_ind1, z_ind2]
        case_dir_map[case]['weights'] = [w1, w2]
    
    return case_dir_map

# In[]:
def prepare_NREL_data (WRF_result_base_loc, vhub):
    NREL_Data = pd.read_csv(os.path.join(WRF_result_base_loc,'NREL_5MW_126_RWT.csv'))
    ws_NREL    = list(NREL_Data['Wind Speed [m/s]'])
    power_NREL = list(NREL_Data['Power [kW]'])

    f = interpolate.interp1d(ws_NREL, power_NREL)
    power_interp = f(vhub)
    
    return NREL_Data, ws_NREL, power_NREL, power_interp

# In[]
def initialize_tabulated_data (power_interp, ind1_for_tab, ind2_for_tab, ind3_for_tab):
    case_tab = ['NREL']
    power1_tab = [power_interp[ind1_for_tab]*1.0e+3]
    power2_tab = [power_interp[ind2_for_tab]*1.0e+3]
    power3_tab = [power_interp[ind3_for_tab]*1.0e+3]
    error1_tab = [0.0]
    error2_tab = [0.0]
    error3_tab = [0.0]

    return case_tab, power1_tab, power2_tab, power3_tab, error1_tab, error2_tab, error3_tab

# In[]:
def read_power_for_a_wind_speed (WRF_result_base_loc, case, outfile, ind_vhub, ws_hub, power_case_ws):
    case_for_ws = path.join(WRF_result_base_loc, case, 'power_curve_{}'.format(ws_hub))
    #print ('case_for_ws: {}'.format(case_for_ws))
    case_nc_file = path.join(case_for_ws, outfile)
    #print ('case_nc_file: {}'.format(case_nc_file))
    case_data = Dataset(case_nc_file, mode='r')
    case_power = case_data.variables['POWER'][:]
    #print ('case_power: ', case_power)
    case_data.close()

    if (len(case_power) <2):
        power_case_ws[ind_vhub] = np.nan
    else:
        power_case_ws[ind_vhub] = case_power[1]
        
    return power_case_ws


# In[]:
def read_power_for_a_case (WRF_result_base_loc, case, vhub, power_interp, outfile, \
                           case_dir_map, case_tab, \
                           power1_tab, power2_tab, power3_tab, \
                           error1_tab, error2_tab, error3_tab, \
                           ind1_for_tab, ind2_for_tab, ind3_for_tab):
    print('Case: {}'.format(case))
    print('Legend: {}'.format(case_dir_map[case]['legend']))
    case_tab.append(case)

    power_case_ws = np.zeros(len(vhub))
    for ind_vhub, ws_hub in enumerate(vhub):
        power_case_ws = read_power_for_a_wind_speed (WRF_result_base_loc, case, outfile, \
                                                     ind_vhub, ws_hub, power_case_ws)
    #print('power_case_ws: ', power_case_ws)
    case_dir_map[case]['power'] = power_case_ws

    error_NREL = (power_case_ws*1e-3/power_interp - 1.0)*100
    case_dir_map[case]['error'] = error_NREL

    print('Power: {}'.format(case_dir_map[case]['power']))
    print('Error w.r.t. NREL: {}'.format(case_dir_map[case]['error']))

    power1_tab.append(float('%.0f'%(power_case_ws[ind1_for_tab])))
    power2_tab.append(float('%.0f'%(power_case_ws[ind2_for_tab])))
    power3_tab.append(float('%.0f'%(power_case_ws[ind3_for_tab])))
    error1_tab.append(float('%5.2f'%(error_NREL[ind1_for_tab])))
    error2_tab.append(float('%5.2f'%(error_NREL[ind2_for_tab])))
    error3_tab.append(float('%5.2f'%(error_NREL[ind3_for_tab])))

    print ('\n')
    
    return case_dir_map, case_tab, power1_tab, power2_tab, power3_tab, error1_tab, error2_tab, error3_tab


# In[]:
def tabulate_data_for_cases (tab_data_file, vhub, case_tab, \
                           power1_tab, power2_tab, power3_tab, \
                           error1_tab, error2_tab, error3_tab, \
                           ind1_for_tab, ind2_for_tab, ind3_for_tab):
    tabulated_data = pd.DataFrame(index = case_tab)
    tabulated_data['power[vhub = {}]'.format(vhub[ind1_for_tab])] = power1_tab
    tabulated_data['error[vhub = {}]'.format(vhub[ind1_for_tab])] = error1_tab
    tabulated_data['power[vhub = {}]'.format(vhub[ind2_for_tab])] = power2_tab
    tabulated_data['error[vhub = {}]'.format(vhub[ind2_for_tab])] = error2_tab
    tabulated_data['power[vhub = {}]'.format(vhub[ind3_for_tab])] = power3_tab
    tabulated_data['error[vhub = {}]'.format(vhub[ind3_for_tab])] = error3_tab
    print('Tabulated Dataframe: \n{}\n'.format(tabulated_data))
    tabulated_data.to_csv(tab_data_file)
    
    
# In[]:
def plot_power_curve_for_cases (NREL_Data, vhub, case_keys, case_dir_map, plt_title, \
                                figFileName, savefig):
    pr = 1.00

    #plt.cla()
    # POWER Curve
    plt.figure()
    plt.plot(NREL_Data['Wind Speed [m/s]'], NREL_Data['Power [kW]']*1e-3/pr, 'k-', label = 'NREL data')
    #plt.plot([0,20],[5,5],'k--',label='Rated power')
    for case in case_keys:
        power_case_ws = case_dir_map[case]['power']
        legend = case_dir_map[case]['legend']
        line_style = case_dir_map[case]['line_style']
        plt.plot(vhub, power_case_ws*1e-6/pr, line_style, label=legend)
    plt.title(plt_title,fontsize=12)
    plt.xlim(2.9,15.1)
    plt.xticks(np.arange(3,15,1.5))
    plt.ylim(0,7.0)
    plt.xlabel(r'Hub-height wind speed [m s$^{-1}$]',fontsize=12)
    plt.ylabel(r'Power (MW)',fontsize=12)
    plt.legend(loc='best')
    plt.tick_params(axis = 'x', labelsize=12)
    plt.tick_params(axis = 'y', labelsize=12)
    
    if savefig:
        plt.savefig(figFileName, bbox_inches='tight')
        
    plt.show()
        
# In[]:
def plot_power_curve_error_for_cases (vhub, case_keys, case_dir_map, plt_title):
    #plt.cla()
    # Error w.r.t. NREL data
    plt.figure()
    for case in case_keys:
        error_case_ws = case_dir_map[case]['error']
        legend = case_dir_map[case]['legend']
        line_style = case_dir_map[case]['line_style']
        plt.plot(vhub, error_case_ws, line_style, label=legend)
    plt.title(plt_title,fontsize=12)
    plt.xlim(2.9,15.1)
    plt.xticks(np.arange(3,15,1.5))
    plt.xlabel(r'Hub-height wind speed [m s$^{-1}$]',fontsize=12)
    plt.ylabel(r'Percent error w.r.t. NREL data',fontsize=12)
    plt.legend(loc='best')
    plt.gca().tick_params(labelsize=10)

    plt.show()
    
    
# In[]:
def read_wrfout_data_for_case_ws (WRF_result_base_loc, case, ws, outfile):
    case_loc = os.path.join(WRF_result_base_loc, case)
    case_ws_loc = os.path.join(case_loc, 'power_curve_{}'.format(ws))
    
    case_ws_outfile = path.join(case_ws_loc, outfile)
    case_ws_wrf_data = xr.open_dataset(case_ws_outfile)
    
    return case_ws_wrf_data

# In[]:
def read_tsout_data_for_case_ws (WRF_result_base_loc, case, ws, tsoutfile):
    case_loc = path.join(WRF_result_base_loc, case)
    case_ws_loc = path.join(case_loc, 'power_curve_{}'.format(ws))
    #print('case/wind_speed loc: {}'.format(case_ws_loc))
    
    case_ws_tsoutfile = path.join(case_ws_loc, tsoutfile)
    case_ws_ts_data = xr.open_dataset(case_ws_tsoutfile)
    
    return case_ws_ts_data

# In[]:
def get_hub_height_slice_data_case_ws (WRF_result_base_loc, case, ws, outfile, tsoutfile, \
                                      time_ind, z_ind1, z_ind2, w1, w2):
    case_ws_wrf_data = read_wrfout_data_for_case_ws (WRF_result_base_loc, case, ws, outfile)
    case_ws_ts_data = read_tsout_data_for_case_ws (WRF_result_base_loc, case, ws, tsoutfile)
    DX, DY = case_ws_wrf_data.DX, case_ws_wrf_data.DY
    #print('DX: {}, DY: {}, DZ: {}'.format(DX, DY, DZ))

    slice_data1_U = extract_slices_from_tsout_file (case_ws_ts_data, 'UTS', z_ind1, time_ind)
    slice_data1_V = extract_slices_from_tsout_file (case_ws_ts_data, 'VTS', z_ind1, time_ind)
    slice_data1_UMAG = compute_u_mag_for_slice (slice_data1_U, slice_data1_V)

    slice_data2_U = extract_slices_from_tsout_file (case_ws_ts_data, 'UTS', z_ind2, time_ind)
    slice_data2_V = extract_slices_from_tsout_file (case_ws_ts_data, 'VTS', z_ind2, time_ind)
    slice_data2_UMAG = compute_u_mag_for_slice (slice_data2_U, slice_data2_V)

    slice_data_wgt_UMAG = w1*slice_data1_UMAG + w2*slice_data2_UMAG
    
    return slice_data_wgt_UMAG, DX, DY
    
    
# In[]:
def get_hub_height_slice_data_in_case_data_map (case_dir_map, WRF_result_base_loc, ws, \
                                                outfile, tsoutfile, time_ind):
    for case in case_dir_map.keys():
        z_ind1 = case_dir_map[case]['z_indices'][0]
        z_ind2 = case_dir_map[case]['z_indices'][1]
        w1 = case_dir_map[case]['weights'][0]
        w2 = case_dir_map[case]['weights'][1]
        #print('z_indices: [{}, {}], weights: [{}, {}]'.format(z_ind1, z_ind2, w1, w2))

        slice_data_wgt_UMAG, DX, DY = \
        get_hub_height_slice_data_case_ws (WRF_result_base_loc, case, ws, outfile, tsoutfile, \
                                      time_ind, z_ind1, z_ind2, w1, w2)
        
        case_dir_map[case]['slice_data'] = slice_data_wgt_UMAG
        case_dir_map[case]['DX'] = DX
        case_dir_map[case]['DY'] = DY
        
    return case_dir_map

# In[]:
def plot_hub_height_slice_data_all_cases_separately (case_dir_map, GAD_param_slice_loc, slice_dir, \
                                                     ws, cont_levels_count, qoi_cont_range, \
                                                     xlim, ylim):
    for case in case_dir_map.keys():
        case_legend = case_dir_map[case]['legend']

        slice_data_wgt_UMAG = case_dir_map[case]['slice_data']
        DX = case_dir_map[case]['DX']
        DY = case_dir_map[case]['DY']

        plot_loc = os.path.join(GAD_param_slice_loc, slice_dir, 'power_curve_{}'.format(ws))
        os.system('mkdir -p %s'%plot_loc)
        image_name = 'slice_{}_{}.png'.format(case, 'UMAG')

        plot_contour_slice(slice_data_wgt_UMAG, DX, DY, plot_loc, image_name, \
                            'UMAG', 'UMAG', 'm/s', case_legend, \
                            cont_levels_count, qoi_cont_range, xlim, ylim)
