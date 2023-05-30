######################################################
#
# Script to view preliminary results
# RSA 6-9-17
#
#####################################################
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

#inputs
savefig = False
figfile = 'power_curve.eps'

#wrf-gad
dirall = './'
outfile = 'wrfout_d01_0001-01-01_00:00:00'
wrfdir =['power_curve_3/', \
         'power_curve_4.5/', \
         'power_curve_6/', \
         'power_curve_7.5/', \
         'power_curve_9/', \
         'power_curve_10.5/', \
         'power_curve_12/', \
         'power_curve_13.5/', \
         'power_curve_15/']
         # 'power_curve_18/']
hub_height = 80.0
vhub = [3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15]

power = np.zeros(len(vhub))
for nd,d in enumerate(wrfdir):
    #get variables from netcdf
    wrfout = dirall + d + outfile
    print wrfout
    data = Dataset(wrfout,mode='r')

    power_d = data.variables['POWER'][:]
    print ('power_d: ', power_d)
    data.close()

    power[nd] = power_d[1]

    # un = data.variables['U'][:,k,:,:]
    # vn = data.variables['V'][:,k,:,:]
    # uc = 0.5 * (un[:,:,1:] + un[:,:,:-1])
    # vc = 0.5 * (vn[:,1:,:] + vn[:,:-1,:])
    # mag = np.sqrt(uc**2 + vc**2)
    # ph = data.variables['PH'][0,:,0,0]
    # phb = data.variables['PHB'][0,:,0,0]
    # zn = ( ph[:-1] + phb[:-1] + ph[1:] + phb[1:] ) / 2 / 9.81
    # print zn[:20]

#actual
#from pdf sent by sonia
print('power: ', power)
va = np.arange(3.0,13.5,0.5)
pa = np.array([0,20,63,116,177,248,331,428,540,667,812,972,1141,1299,1448,1561,1633,1661,1677,1678,1680])/1000.0 #MW

#plot
pr = 1.68

plt.figure()
plt.plot(va,pa/pr,'k-',label='Actual')
plt.plot(vhub,power*1e-6/pr,'m.-',label='WRF-GAD model')
plt.plot([0,20],[1,1],'k--',label='Rated power')
plt.xlim(0,18)
plt.xticks(np.arange(0,21,3))
plt.ylim(0,1.2)
plt.xlabel(r'Hub-height Wind Speed [m s$^{-1}$]',fontsize=12)
plt.ylabel(r'$P/P_R$',fontsize=12)
plt.legend(loc='lower right')
plt.gca().tick_params(labelsize=10)

if savefig:
    plt.savefig(figfile,dpi=100)
    plt.show()
else:
    plt.show()
