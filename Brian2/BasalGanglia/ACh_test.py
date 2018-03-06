# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 18:07:54 2018

@author: Elijah
"""

from __future__ import division
from brian2 import *
import numpy as np
prefs.codegen.target = 'cython'

start_scope()     

n = 20

# TC Iz
c = -63
b= 0.25 
a = 0.02
d = 0.05

# plasticity
traceTau = traceTauPost = 1*second  # this will change integration window 
traceConstant = 0.01 # this will change maximum weight change
tracePostConstant = -traceConstant*traceTau/traceTauPost*1.05

wmax = 10

# Equations
eqs = '''
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = a*(b*v-u)/ms : 1
I : 1
''' 

reset='''
v=c
u += d
'''

# ACh plasticity
ACh_stdp = '''
             w : 1
             dtraceCon/dt = -traceCon/traceTau : 1 (event-driven)
             dtraceConPost/dt = -traceConPost/traceTauPost : 1 (event-driven)
             '''   
onpre_ACh = '''
             v_post += w
             traceCon += traceConstant
             w = clip(w-traceConPost, 0, wmax)
             '''                                              
onpost_ACh = '''
             traceConPost += tracePostConstant 
             w = clip(w-traceCon, 0, wmax)
             '''

## Neurons             
ACh = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset,method='euler')

ACh.v = c
ACh.u = b*c
ACh.I = 0


Cortex = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset,method='euler')

Cortex.v = c
Cortex.u = b*c
Cortex.I = 0


MSN = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset,method='euler')

MSN.v = c
MSN.u = b*c
MSN.I = 0

## Connections

ACh_striatum_prob = 40/n
ACh_striatum = Synapses(ACh, MSN, ACh_stdp, on_pre=onpre_ACh, on_post=onpost_ACh) #on_pre='v=+10')
ACh_striatum.connect(j='k for k in range(i-n, i+n) if rand()<ACh_striatum_prob', skip_if_invalid=True)
ACh_striatum.delay = 10*ms 
ACh_striatum.w = 5 

Cortex_MSN = Synapses(Cortex, MSN, ACh_stdp, on_pre=onpre_ACh, on_post=onpost_ACh) #on_pre='v=+10')
Cortex_MSN.connect(j='k for k in range(i-n, i+n) if rand()<ACh_striatum_prob', skip_if_invalid=True)
Cortex_MSN.delay = 10*ms 
Cortex_MSN.w = 5 


def ACh_plasticity(t,window,SpikeMon,SpikeMon2,SpikeMon3,SynapseMon,SynapseMon2):          
    AChind = np.where(SpikeMon2.t > (t - window)) 
    if len(AChind[0]) > 5: # was DA released?
       CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
       PresynapticInd = []
       for i in range(0,len(np.unique(CortexIndex))):
           pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
           PresynapticInd.extend(pre[0])
       MSN_High = SpikeMon3.i    
       act_MSN = np.where(SpikeMon3 > MSN_High) # find MSNs that are in a "high" state
       high_synaptic_ind = np.concatenate([np.where(SynapseMon.j == x)[0] for x in act_MSN[0]]) # synaptic indicies projecting to high-state MSNs
       s2 = set(high_synaptic_ind) # set object for high_synpatic ind
       strengthen_synapse = [val for val in PresynapticInd if val in s2] # strengthen synapses that have MSNs in up state and cortical/DA input
       SynapseMon.w[strengthen_synapse] -=  10*(SynapseMon.traceCon[strengthen_synapse] * mean(SynapseMon2.traceCon))     
       SynapseMon.w = clip(SynapseMon.w, 0, wmax)
       
 
rew_win = 5*ms  
cortex_spikes = SpikeMonitor(Cortex)    
AChSpikes = SpikeMonitor(ACh)
MSNspikes = SpikeMonitor(MSN)

@network_operation(dt=rew_win)
def calculate_LTP(t):
    ACh_plasticity(t,rew_win, cortex_spikes, AChSpikes, MSNspikes, Cortex_MSN, ACh_striatum)    


# Run the sequence 
ACh.I = 0
Cortex.I = 10
MSN.I = 10
run(1000*ms,report='text')
ACh.I = -100
run(1000*ms, report='text')
ACh.I = 0
run(1000*ms, report='text')          
          
figure()
plot(AChSpikes.t, AChSpikes.i)
xlabel('Time(ms)')
ylabel('Firing Rate')
title('ACh Firing Rates')
         
       