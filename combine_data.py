#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
def combine_power_for_intervals(pickled_data_combined, pickled_data_for_interval):
    power_avg_key = 'power_avg'
    power_inst_key = 'power_inst'
    
    # Append the averaged power
    if power_avg_key not in pickled_data_combined:
        pickled_data_combined[power_avg_key] = pickled_data_for_interval[power_avg_key]
    else:
        pickled_data_combined[power_avg_key] = np.vstack((pickled_data_combined[power_avg_key], pickled_data_for_interval[power_avg_key]))
    
    # Append the instantaneous power
    if power_inst_key not in pickled_data_combined:
        pickled_data_combined[power_inst_key] = pickled_data_for_interval[power_inst_key]
    else:
        pickled_data_combined[power_inst_key] = np.vstack((pickled_data_combined[power_inst_key], pickled_data_for_interval[power_inst_key]))
    
    # Return the combined data
    return pickled_data_combined