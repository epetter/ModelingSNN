# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 18:58:29 2017

@author: Elijah
"""

import numpy as np

def spike_covariance(**kwargs):
    # this function will calculat covariance between n firing rates
        # example
    #spike_cov, firing_rates = spike_covariance(f1=SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,
    #                             f2=D1pop.smooth_rate(window='gaussian',width=binSize)/Hz)
    return np.cov(np.array([x for key, x in kwargs.iteritems()])), np.array([x for key, x in kwargs.iteritems()])
    
a,b = spike_covariance(SNrL=SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,
                 SNrNL=SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,
                 SNrWL=SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,
                 STN=STNpop.smooth_rate(window='gaussian',width=binSize)/Hz,
                 D1pop=D1pop.smooth_rate(window='gaussian',width=binSize)/Hz)
                 
plt.plot(np.transpose(b)) 
print a                
                 
    

spike_covariance(SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,STNpop.smooth_rate(window='gaussian',width=binSize)/Hz    

spike_covariance(SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,D1pop.smooth_rate(window='gaussian',width=binSize)/Hz