# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 18:58:29 2017

@author: Elijah
"""

import numpy as np

def spike_covariance(firing_rate1, firing_rate2):
    return np.cov(np.array([firing_rate1,firing_rate2]))
    
    