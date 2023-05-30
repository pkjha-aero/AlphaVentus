######################################################
#
# Script to view preliminary results
# RSA 6-9-17
#
#####################################################
import sys
import os.path as path

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

#inputs
savefig = False
figfile = 'power_curve.eps'

#wrf-gad
dirall_case1 = sys.argv[1]
dirall_case2 = sys.argv[2]

legend1 = dirall_case1
legend2 = dirall_case2

outfile = 'wrfout_d01_0001-01-01_00:00:00'
wrfdir =['power_curve_3', \
         'power_curve_4.5', \
         'power_curve_6', \
         'power_curve_7.5', \
         'power_curve_9', \
         'power_curve_10.5', \
         'power_curve_12', \
         'power_curve_13.5', \
         'power_curve_15', \
         'power_curve_18']
hub_height = 80.0
#vhub = [4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 18]
vhub = [3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 18]

# NREL data
NREL_Data = pd.read_csv('NREL_5MW_126_RWT.csv')
ws_NREL    = list(NREL_Data['Wind Speed [m/s]'])
power_NREL = list(NREL_Data['Power [kW]'])

f = interpolate.interp1d(ws_NREL, power_NREL)
vhub = [3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 18]
power_interp = f(vhub)

#plt.plot(NREL_Data['Wind Speed [m/s]'], NREL_Data['Power [kW]'], label = 'NREL data')
#plt.plot(vhub, power_interp, '.-', label = 'NREL interpolated')


# Extract Sim Data

power_case1 = np.zeros(len(vhub))
power_case2 = np.zeros(len(vhub))
for nd,d in enumerate(wrfdir):
    #get variables from netcdf
    wrfout_case1 = path.join(dirall_case1, d, outfile)
    wrfout_case2 = path.join(dirall_case2, d, outfile)
    print ('wrfout_case1', wrfout_case1)
    print ('wrfout_case2', wrfout_case2)
    data_case1 = Dataset(wrfout_case1,mode='r')
    data_case2 = Dataset(wrfout_case2,mode='r')

    power_d_case1 = data_case1.variables['POWER'][:]
    power_d_case2 = data_case2.variables['POWER'][:]
    print ('power_d_case1: ', power_d_case1)
    print ('power_d_case2: ', power_d_case2)
    data_case1.close()
    data_case2.close()
    
    if (len(power_d_case1) < 2):
        power_case1[nd] = np.nan
    else:
        power_case1[nd] = power_d_case1[1]
        
    if (len(power_d_case2) < 2):
        power_case2[nd] = np.nan
    else:
        power_case2[nd] = power_d_case2[1]

#actual
#from pdf sent by sonia
print('power_case1: ', power_case1)
print('power_case2: ', power_case2)

percent_error = (power_case2/power_case1 -1.0)*100  # Error in case 2 w.r.t. case 1
percent_error1 = (power_case1*1e-3/power_interp - 1.0)*100  # Error in case 1 w.r.t. interpolated NREL data
percent_error2 = (power_case2*1e-3/power_interp - 1.0)*100  # Error in case 2 w.r.t. interpolated NREL data

#va = np.arange(3.0,13.5,0.5)
#pa = np.array([0,20,63,116,177,248,331,428,540,667,812,972,1141,1299,1448,1561,1633,1661,1677,1678,1680])/1000.0 #MW

#plot
pr = 1.00

# POWER Curve
plt.figure(1)
#plt.plot(va,pa/pr,'k-',label='Actual')
plt.plot(vhub,power_case1*1e-6/pr,'b--',label=legend1)
plt.plot(vhub,power_case2*1e-6/pr,'m.-',label=legend2)
plt.plot(NREL_Data['Wind Speed [m/s]'], NREL_Data['Power [kW]']*1e-3/pr, 'k-', label = 'NREL data')
plt.plot([0,20],[5,5],'k--',label='Rated power')
plt.xlim(0,18)
plt.xticks(np.arange(0,21,3))
plt.ylim(0,7.0)
plt.xlabel(r'Hub-height Wind Speed [m s$^{-1}$]',fontsize=12)
plt.ylabel(r'P (MW)',fontsize=12)
plt.legend(loc='lower right')
plt.gca().tick_params(labelsize=10)

# Error w.r.t. NREL data
plt.figure(2)
plt.plot(vhub,percent_error1,'b--',label='%s w.r.t. %s'%(legend1, 'NREL Data'))
plt.plot(vhub,percent_error2,'m.-',label='%s w.r.t. %s'%(legend2, 'NREL Data'))
plt.xlim(0,18)
plt.xticks(np.arange(0,21,3))
plt.xlabel(r'Hub-height Wind Speed [m s$^{-1}$]',fontsize=12)
plt.ylabel(r'Percent Error',fontsize=12)
plt.legend(loc='lower right')
plt.gca().tick_params(labelsize=10)

plt.figure(3)
plt.plot(vhub,percent_error,'b--',label='%s w.r.t. %s'%(legend2, legend1))
plt.xlim(0,18)
plt.xticks(np.arange(0,21,3))
plt.xlabel(r'Hub-height Wind Speed [m s$^{-1}$]',fontsize=12)
plt.ylabel(r'Percent Error',fontsize=12)
plt.legend(loc='lower right')
plt.gca().tick_params(labelsize=10)

if savefig:
    plt.savefig(figfile,dpi=100)
    plt.show()
else:
    plt.show()
