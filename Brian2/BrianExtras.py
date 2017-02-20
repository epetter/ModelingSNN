# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 11:10:44 2017

@author: elijah
"""

from brian2 import *
import numpy as np
#%%
# Extra Brain Functions 

def visualise_connectivity(S):
    Ns = len(S.source)
    Nt = len(S.target)
    figure(figsize=(10, 4))
    subplot(121)
    plot(zeros(Ns), arange(Ns), 'ok', ms=10)
    plot(ones(Nt), arange(Nt), 'ok', ms=10)
    for i, j in zip(S.i, S.j):
        plot([0, 1], [i, j], '-k')
    xticks([0, 1], ['Source', 'Target'])
    ylabel('Neuron index')
    xlim(-0.1, 1.1)
    ylim(-1, max(Ns, Nt))
    subplot(122)
    plot(S.i, S.j, 'ok')
    xlim(-1, Ns)
    ylim(-1, Nt)
    xlabel('Source neuron index')
    ylabel('Target neuron index')
    
def frange(start, stop, step):
 i = start
 while i < stop:
     yield i
     i += step
     
def calculate_FR(SpikeMon,binSize=100*ms,timeWin=learn_duration):
    allBin = []
    for i in frange(360,int(ceil(timeWin/ms)),binSize/ms):
         binz = len(np.where((SpikeMon.t > i*ms)*(SpikeMon.t < (i*ms+binSize)))[0])/binSize
         allBin.append(binz)
    FR = np.mean(allBin)
    return FR, allBin       
       