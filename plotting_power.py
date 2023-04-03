#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 19:14:28 2023

@author: jha3
"""

# In[]
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
from datetime import date, datetime, timedelta, time
import pickle

# In[]
def plot_power_inst (pickle_file_name, plot_loc, case_name, dt, ylim=None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
        
    plt.figure()
    time_sec = np.arange(len(pickled_data_read['power_inst']))*dt
    time_min = time_sec/60.0
    plt.plot(time_min, pickled_data_read['power_inst']/1.0e6)

    plt.xlabel('Time [Min]', fontsize=14)
    plt.ylabel('Power [MW]', fontsize=14)
    if ylim:
        plt.ylim([ylim[0],ylim[1]])
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title('Time series of power', fontsize=14)

    if ylim:
        filename = 'Power_TS_{}_Bounded.png'.format(case_name)
    else:
        filename = 'Power_TS_{}_Unbounded.png'.format(case_name)
    filedir = os.path.join(plot_loc, 'Instantaneous', 'Power')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()
     

# In[]
def plot_power_pdf (pickle_file_name, plot_loc, case_name, num_bins = 20, xlim = None, ylim = None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
        
    plt.figure()
    plt.bar(pickled_data_read['bin_centers'], pickled_data_read['hist'], color = 'lime')
    plt.plot(pickled_data_read['bin_centers'], pickled_data_read['hist'], color = 'k')

    plt.ylabel('PDF of power[1/MW]', fontsize=14)
    plt.xlabel('Power [MW]', fontsize=14)
    if xlim:
        plt.xlim([xlim[0],xlim[1]])
    if ylim:
        plt.ylim([ylim[0],ylim[1]])
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title('Probability density function (PDF) of power', fontsize=14)

    if xlim or ylim:
        filename = 'Power_PDF_{}_Bounded.png'.format(case_name)
    else:
        filename = 'Power_PDF_{}_Unbounded.png'.format(case_name)
    filedir = os.path.join(plot_loc, 'Instantaneous', 'Power')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()
    
# In[]
def plot_power_pdf_all_cases (pickle_file_name, plot_loc, case_colors, case_names, group_identifier, xlim = None, ylim = None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
     
    plt.figure()
    
    colors = ['r', 'b', 'g', 'm']
    for case_name  in case_names:
        plt.plot(pickled_data_read[case_name]['bin_centers'], pickled_data_read[case_name]['hist'], color = case_colors[case_name], label=case_name)
    
    plt.ylabel('PDF of power [1/MW]', fontsize=14)
    plt.xlabel('Power [MW]', fontsize=14)
    if xlim:
        plt.xlim([xlim[0],xlim[1]])
    if ylim:
        plt.ylim([ylim[0],ylim[1]])
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title('Probability density function (PDF) of power', fontsize=14)
    plt.legend(fontsize=14, ncol = 1)
    #plt.legend(bbox_to_anchor=(1.0, 0.7), fontsize=14, ncol = 1)
    
    if xlim or ylim:
        filename = 'Power_PDF_{}_Bounded.png'.format(group_identifier)
    else:
        filename = 'Power_PDF_{}_Unbounded.png'.format(group_identifier)
    filedir = os.path.join(plot_loc, 'Instantaneous', 'Power')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()
    
    
# In[]
def plot_power_inst_pdf (pickle_file_name, plot_loc, case_name, dt, num_bins = 20, bounds_for_power = None, bounds_for_pdf = None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
    
    fig, ax = plt.subplots(1, 2, sharey = True)
    
    # In[]: Time Series
    time_sec = np.arange(len(pickled_data_read['power_inst']))*dt
    time_min = time_sec/60.0
    ax[0].plot(time_min, pickled_data_read['power_inst']/1.0e6)

    ax[0].set_xlabel('Time [Min]', fontsize=14)
    ax[0].set_ylabel('Power [MW]', fontsize=14)
    if bounds_for_power:
        ax[0].set_ylim([bounds_for_power[0],bounds_for_power[1]])
    ax[0].tick_params(axis='x', labelsize=14)
    ax[0].tick_params(axis='y', labelsize=14)
    
    # In[]: PDF
    ax[1].barh(pickled_data_read['bin_centers'], pickled_data_read['hist'], color = 'lime')
    ax[1].plot(pickled_data_read['hist'], pickled_data_read['bin_centers'], color = 'k')

    ax[1].set_xlabel('PDF of power[1/MW]', fontsize=14)
    if bounds_for_power:
        ax[1].set_ylim([bounds_for_power[0],bounds_for_power[1]])
    if bounds_for_pdf:
        ax[1].set_xlim([bounds_for_pdf[0],bounds_for_pdf[1]])
    ax[1].tick_params(axis='x', labelsize=14)
    ax[1].tick_params(axis='y', labelsize=14)

    fig.suptitle('Instantaneous power and its PDF', fontsize=14)
    
    # In[]: Filename
    filename = 'Power_TS_PDF_{}_Unbounded.png'.format(case_name)
    if bounds_for_power and not bounds_for_pdf:
        filename = 'Power_TS_PDF_{}_Bounded_Power.png'.format(case_name)
    if bounds_for_pdf and not bounds_for_power:
        filename = 'Power_TS_PDF_{}_Bounded_PDF.png'.format(case_name)
    if bounds_for_power and bounds_for_pdf:
        filename = 'Power_TS_PDF_{}_Bounded_Power_PDF.png'.format(case_name)
        
    filedir = os.path.join(plot_loc, 'Instantaneous', 'Power')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()


# In[]
def plot_power_spectra_combined (pickle_file_name, plot_loc, case_name, xlim = None, ylim = None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
        
    plt.figure()
    plt.loglog(pickled_data_read['freq'], pickled_data_read['psd'])
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('PSD [MW^2/Hz]')
    if xlim:
        plt.xlim([xlim[0],xlim[1]])
    if ylim:
        plt.ylim([ylim[0],ylim[1]])
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title('Power spectral density of turbine power', fontsize=14)

    if xlim or ylim:
        filename = 'Power_PSD_{}_Bounded.png'.format(case_name)
    else:
        filename = 'Power_PSD_{}_Unbounded.png'.format(case_name)
    filedir = os.path.join(plot_loc, 'Instantaneous', 'Power')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()


# In[]
def plot_power_spectra_all_cases (pickle_file_name, plot_loc, case_colors, case_names, group_identifier, xlim = None, ylim = None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
     
    plt.figure()
    
    for case_name  in case_names:
        plt.loglog(pickled_data_read[case_name]['freq'], pickled_data_read[case_name]['psd'], color = case_colors[case_name], label=case_name)
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('PSD [MW^2/Hz]')
    if xlim:
        plt.xlim([xlim[0],xlim[1]])
    if ylim:
        plt.ylim([ylim[0],ylim[1]])
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title('Power spectral density of turbine power', fontsize=14)
    plt.legend(fontsize=14, ncol = 1)
    #plt.legend(bbox_to_anchor=(1.0, 0.7), fontsize=14, ncol = 1)
    
    if xlim or ylim:
        filename = 'Power_PSD_{}_Bounded.png'.format(group_identifier)
    else:
        filename = 'Power_PSD_{}_Unbounded.png'.format(group_identifier)
    filedir = os.path.join(plot_loc, 'Instantaneous', 'Power')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()

# In[]
def plot_power_avg (pickle_file_name, plot_loc, case_name, ylim=None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
        
    power_avg = pickled_data_read['power_avg']
    intervals = np.array(range(len(power_avg))) + 1
    plt.figure()
    plt.bar(intervals, power_avg.flatten()/1.0e6, width = 0.4)

    plt.xlabel('Interval', fontsize=14)
    plt.ylabel('Power [MW]', fontsize=14)
    if ylim:
        plt.ylim([ylim[0],ylim[1]])
    plt.xticks(intervals)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title('Average power over 10-min', fontsize=14)

    if ylim:
        filename = 'Power_TimeAvg_{}_Bounded.png'.format(case_name)
    else:
        filename = 'Power_TimeAvg_{}_Unbounded.png'.format(case_name)
    filedir = os.path.join(plot_loc, 'TimeAvg', 'Power_Avg')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()


# In[]
def plot_power_stdev_combined (pickle_file_name, plot_loc, case_name, ylim=None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
        
    power_stdev = pickled_data_read['power_stdev']
    intervals = np.array(range(len(power_stdev))) + 1
    plt.figure()
    plt.bar(intervals, power_stdev.flatten()/1.0e6, width = 0.4)

    plt.xlabel('Interval', fontsize=14)
    plt.ylabel('Power [MW]', fontsize=14)
    if ylim:
        plt.ylim([ylim[0],ylim[1]])
    plt.xticks(intervals)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title('Std. deviation in power over 10-min', fontsize=14)

    if ylim:
        filename = 'Power_StDev_{}_Bounded.png'.format(case_name)
    else:
        filename = 'Power_StDev_{}_Unbounded.png'.format(case_name)
    filedir = os.path.join(plot_loc, 'TimeAvg', 'Power_Avg')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()
    
# In[]
def plot_power_avg_all_cases (pickle_file_name, plot_loc, ylim=None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
        
    case_names = list(pickled_data_read.keys())
    power_avg = pickled_data_read[case_names[0]]['power_avg']
    intervals = np.array(range(len(power_avg))) + 1
    
    power_avg_all_cases = {}
    
    for case_name in case_names:
        power_avg_all_cases[case_name] = pickled_data_read[case_name]['power_avg'].flatten()/1.0e6
    
    df = pd.DataFrame(power_avg_all_cases, index = intervals)
    
    plt.figure()
    ax = df.plot.bar(rot=0)

    ax.set_xlabel('Interval', fontsize=14)
    ax.set_ylabel('10-min averaged power [MW]', fontsize=14)
    if ylim:
        ax.set_ylim([ylim[0],ylim[1]])

    #ax.set_xticks(intervals)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    #plt.title('Average power over 10-min', fontsize=14)
    ax.legend(bbox_to_anchor=(1.0, 0.7), fontsize=14, ncol = 1)
    
    if ylim:
        filename = 'Power_TimeAvg_{}_Bounded.png'.format(case_name)
    else:
        filename = 'Power_TimeAvg_{}_Unbounded.png'.format(case_name)
    filedir = os.path.join(plot_loc, 'TimeAvg', 'Power_Avg')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()
    
# In[]
def plot_power_stdev_all_cases (pickle_file_name, plot_loc, ylim=None):
    # In[] Read the pickle data
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
        
    case_names = list(pickled_data_read.keys())
    power_stdev = pickled_data_read[case_names[0]]['power_stdev']
    intervals = np.array(range(len(power_stdev))) + 1
    
    power_stdev_all_cases = {}
    
    for case_name in case_names:
        power_stdev_all_cases[case_name] = pickled_data_read[case_name]['power_stdev'].flatten()/1.0e6
    
    df = pd.DataFrame(power_stdev_all_cases, index = intervals)
    
    plt.figure()
    ax = df.plot.bar(rot=0)

    ax.set_xlabel('Interval', fontsize=14)
    ax.set_ylabel('Std. dev. in power for 10-min interval [MW]', fontsize=14)
    if ylim:
        ax.set_ylim([ylim[0],ylim[1]])

    #ax.set_xticks(intervals)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    #plt.title('Average power over 10-min', fontsize=14)
    ax.legend(bbox_to_anchor=(1.0, 0.7), fontsize=14, ncol = 1)
    
    if ylim:
        filename = 'Power_StDev_{}_Bounded.png'.format(case_name)
    else:
        filename = 'Power_StDev_{}_Unbounded.png'.format(case_name)
    filedir = os.path.join(plot_loc, 'TimeAvg', 'Power_StDev')
    os.system('mkdir -p %s'%filedir)
    plt.savefig(os.path.join(filedir, filename), bbox_inches='tight')
    plt.close()