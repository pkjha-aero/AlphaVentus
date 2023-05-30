######################################################
#
# Script to plot power curve
# PKJ: 2023/01/05
#
#####################################################
import sys
import os
import os.path as path

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#import openpyxl

# Paths
WRF_result_base_loc ='/p/lustre2/jha3/fromWill/TurbTest/NREL_5MW'
GAD_param_output_loc = os.path.join(WRF_result_base_loc, 'GAD_Param_Study_Output')
GAD_param_power_curve_loc = os.path.join(GAD_param_output_loc, 'PowerCurve')
os.system('mkdir -p %s'%GAD_param_power_curve_loc)
#inputs
savefig = True

# Flags
induction_effect = True
grid_effect = False
epsilon_effect = False


#wrf-gad
# Effect of Induction Type
if induction_effect:
	figFileName = os.path.join(GAD_param_power_curve_loc, 'power_induction.png')
	tab_data_file = os.path.join(GAD_param_power_curve_loc, 'power_induction.csv')
	plt_title   = '$Grid: \Delta x = \Delta y = 4 m, \Delta z = 4 m; Gaussian: \epsilon/\Delta_{grid} = 1.00 $'
	case_dir_map = {'Glauert_Tip_dx_04_dz_04_an_fixed_at_iter_eps_1.00': {'legend':'Glauert, a$_{n}$ = 0.07, a$_{t}$ = iterated', 'line_style': 'b'},
		        'Glauert_Tip_dx_04_dz_04_an_iter_at_iter_eps_1.00': {'legend':'Glauert, a$_{n}$ = iterated, a$_{t}$ = iterated', 'line_style': 'g'},
			'Shen_Tip_dx_04_dz_04_an_fixed_at_iter_eps_1.00': {'legend':'Shen, a$_{n}$ = 0.07, a$_{t}$ = iterated', 'line_style': 'r--'},
			'Shen_Tip_dx_04_dz_04_an_iter_at_iter_eps_1.00': {'legend':'Shen, a$_{n}$ = iterated, a$_{t}$ = iterated', 'line_style': 'm--'},
			'Tip_dx_04_dz_04_an_0.0_at_0.0_eps_1.00': {'legend':'a$_{n}$ = 0.0, a$_{t}$ = 0.0', 'line_style': 'c'}}

if grid_effect:
	figFileName = os.path.join(GAD_param_power_curve_loc, 'power_grid_effect.png')
	tab_data_file = os.path.join(GAD_param_power_curve_loc, 'power_grid_effect.csv')
	plt_title   = '$ Induction: Glauert, a_{n} = 0.07, a_{t} = iterated; Gaussian: \epsilon/\Delta_{grid} = 1.00 $'
	case_dir_map = {'Glauert_Tip_dx_04_dz_04_an_fixed_at_iter_eps_1.00': {'legend':'$\Delta x = \Delta y = 4 m, \Delta z = 4 m $', 'line_style': 'b'},
			'Glauert_Tip_dx_08_dz_08_an_fixed_at_iter_eps_1.00': {'legend':'$\Delta x = \Delta y = 8 m, \Delta z = 8 m $', 'line_style': 'g'},
			'Glauert_Tip_dx_16_dz_04_an_fixed_at_iter_eps_1.00': {'legend':'$\Delta x = \Delta y = 16 m, \Delta z = 4 m $', 'line_style': 'r--'},
			'Glauert_Tip_dx_16_dz_16_an_fixed_at_iter_eps_1.00': {'legend':'$\Delta x = \Delta y = 16 m, \Delta z = 16 m $', 'line_style': 'm--'}}


if epsilon_effect:
	figFileName = os.path.join(GAD_param_power_curve_loc, 'power_gaussian_effect.png')
	tab_data_file = os.path.join(GAD_param_power_curve_loc, 'power_gaussian_effect.csv')
	plt_title   = '$Grid: \Delta x = \Delta y = 4 m, \Delta z = 4 ; Induction: Glauert, a_{n} = 0.07, a_{t} = iterated $'
	case_dir_map = {'Glauert_Tip_dx_04_dz_04_an_fixed_at_iter_eps_0.707': {'legend':'$ \epsilon/\Delta_{grid} = 1/\sqrt{2} $', 'line_style': 'b'},
			'Glauert_Tip_dx_04_dz_04_an_fixed_at_iter_eps_1.00': {'legend':'$ \epsilon/\Delta_{grid} = 1.00$', 'line_style': 'g'},
			'Glauert_Tip_dx_04_dz_04_an_fixed_at_iter_eps_1.4142': {'legend':'$ \epsilon/\Delta_{grid} = \sqrt{2} $', 'line_style': 'r--'},
			'Glauert_Tip_dx_04_dz_04_an_fixed_at_iter_eps_2.00': {'legend':'$ \epsilon/\Delta_{grid} = 2.00 $', 'line_style': 'm--'}}


case_keys = list(case_dir_map.keys())
case_keys.sort()
print ('Cases: \n {}'.format(case_keys))

# Common stuff
vhub = [3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15]
#vhub = [3, 4.5, 6, 7.5, 9, 10.5, 13.5, 15, 18]
ind1_for_tab = 3 # To be used in the paper
ind2_for_tab = 4 # To be used in the paper
ind3_for_tab = 5 # To be used in the paper
outfile = 'wrfout_d01_0001-01-01_00:00:00'

# NREL data
NREL_Data = pd.read_csv(os.path.join(WRF_result_base_loc,'NREL_5MW_126_RWT.csv'))
ws_NREL    = list(NREL_Data['Wind Speed [m/s]'])
power_NREL = list(NREL_Data['Power [kW]'])

f = interpolate.interp1d(ws_NREL, power_NREL)
power_interp = f(vhub)

# Read the case data
#power_sim = []

case_tab = ['NREL']
power1_tab = [power_interp[ind1_for_tab]*1.0e+3]
power2_tab = [power_interp[ind2_for_tab]*1.0e+3]
power3_tab = [power_interp[ind3_for_tab]*1.0e+3]
error1_tab = [0.0]
error2_tab = [0.0]
error3_tab = [0.0]

print ('Cases under consideration : \n')
for case in case_keys:
    print('Case: {}'.format(case))
    print('Legend: {}'.format(case_dir_map[case]['legend']))
    case_tab.append(case)

    power_case_ws = np.zeros(len(vhub))
    for ind_vhub, ws_hub in enumerate(vhub):
        case_for_ws = path.join(WRF_result_base_loc, case, 'power_curve_{}'.format(ws_hub))
        #print ('case_for_ws: {}'.format(case_for_ws))
        case_nc_file = path.join(case_for_ws, outfile)
        print ('case_nc_file: {}'.format(case_nc_file))
        case_data = Dataset(case_nc_file, mode='r')
        case_power = case_data.variables['POWER'][:]
	#print ('case_power: ', case_power)
        case_data.close()
        
        if (len(case_power) <2):
            power_case_ws[ind_vhub] = np.nan
        else:
            power_case_ws[ind_vhub] = case_power[1]

    print('power_case_ws: ', power_case_ws)
    case_dir_map[case]['power'] = power_case_ws
    #power_sim.append(power_case_ws)

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

tabulated_data = pd.DataFrame(index = case_tab)
tabulated_data['power[vhub = {}]'.format(vhub[ind1_for_tab])] = power1_tab
tabulated_data['error[vhub = {}]'.format(vhub[ind1_for_tab])] = error1_tab
tabulated_data['power[vhub = {}]'.format(vhub[ind2_for_tab])] = power2_tab
tabulated_data['error[vhub = {}]'.format(vhub[ind2_for_tab])] = error2_tab
tabulated_data['power[vhub = {}]'.format(vhub[ind3_for_tab])] = power3_tab
tabulated_data['error[vhub = {}]'.format(vhub[ind3_for_tab])] = error3_tab
print('Tabulated Dataframe: \n{}\n'.format(tabulated_data))
tabulated_data.to_csv(tab_data_file)
  
#plot
pr = 1.00

plt.cla()
# POWER Curve
plt.figure(1)
plt.plot(NREL_Data['Wind Speed [m/s]'], NREL_Data['Power [kW]']*1e-3/pr, 'k-', label = 'NREL data')
#plt.plot([0,20],[5,5],'k--',label='Rated power')
for case in case_keys:
    power_case_ws = case_dir_map[case]['power']
    legend = case_dir_map[case]['legend']
    line_style = case_dir_map[case]['line_style']
    plt.plot(vhub, power_case_ws*1e-6/pr, line_style, label=legend)
plt.title(plt_title,fontsize=12)
plt.xlim(2,16)
plt.xticks(np.arange(2,16,2))
plt.ylim(0,7.0)
plt.xlabel(r'Hub-height wind speed [m s$^{-1}$]',fontsize=12)
plt.ylabel(r'Power (MW)',fontsize=12)
plt.legend(loc='best')
plt.gca().tick_params(labelsize=12)

if savefig:
    plt.savefig(figFileName,dpi=300)

# Error w.r.t. NREL data
plt.figure(2)

for case in case_keys:
    error_case_ws = case_dir_map[case]['error']
    legend = case_dir_map[case]['legend']
    line_style = case_dir_map[case]['line_style']
    plt.plot(vhub, error_case_ws, line_style, label=legend)
plt.title(plt_title,fontsize=12)
plt.xlim(2,16)
plt.xticks(np.arange(2,16,2))
plt.xlabel(r'Hub-height wind speed [m s$^{-1}$]',fontsize=12)
plt.ylabel(r'Percent error w.r.t. NREL data',fontsize=12)
plt.legend(loc='best')
plt.gca().tick_params(labelsize=10)

plt.show()
