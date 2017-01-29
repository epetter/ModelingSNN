# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 20:59:01 2016

@author: Elijah
"""
from brian2 import*

Npost = 
postTrace = StateMonitor(Npost,'v',record=True)
S = synapses
Smon = StateMonitor(S,'w',record=True)

@network_operation(dt=500*ms)
def update_w:
    maxt = postTrace[max(postTrace)].t # get time? prob not right
    filterz = Smon.t - maxT # subtract time
    Smon.w = filterz*operationz # update weights 
    
    
    
    