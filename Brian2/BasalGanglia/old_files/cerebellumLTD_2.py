# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:54:48 2016

@author: Elijah
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 15:47:32 2016

@author: Elijah
"""

#import numpy as np
from brian2 import *


start_scope()

N_GR = 10
N_PKJ = 10
N_IO = 10

# variables
# thalamo-cortical Iz
c = -66*mV
b= 0.2/ms
a = 0.02/ms
d = 8*mV/ms

c2 = -66*mV
b2= 0.2/ms
a2 = 0.1/ms
d2 = 2*mV/ms


neuron_spacing1 = 50*umetre
width = N_IO/4.0*neuron_spacing1

window = 15*ms

# eqs
eqs =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a*(b*v-u) : volt/second
I : amp
x : metre
active : 1 
'''

reset='''
v=c
u += d
'''

# eqs2
eqs2 =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a2*(b2*v-u) : volt/second
I : amp
x : metre
active : 1 
'''

reset='''
v=c2
u += d2
'''

# Neuron groups
GR = NeuronGroup(N_GR,eqs,threshold='v>30*mV',reset=reset,method='euler') # Granule cells -- their axons are PFs
PKJ = NeuronGroup(N_PKJ,eqs,threshold='v>30*mV',reset=reset,method='euler') #
#PKJ.x = 'i*neuron_spacing1'
IO = NeuronGroup(N_IO,eqs2,threshold='v>40*mV',reset=reset,method='euler') # Inferior Olive neurons -- their axons are CFs 
IO.x = 'i*neuron_spacing1'

# Synapses
CF_PKJ = Synapses(IO, PKJ,'w:1',on_pre='v+=w*mV')
CF_PKJ.connect(condition='i!=j',p=.5)
CF_PKJ.w = '(10*exp(-(x_pre-x_post)**2/(2*width**2)))'

PF_PKJ = Synapses(GR,PKJ,'w:1',on_pre='v+=w*mV')
PF_PKJ.connect(condition='i!=j',p=0.5)# GR -> PKJ randomly, sparsely
PF_PKJ.w = 15 

# need to make the index the event index
GRmon = SpikeMonitor(GR)
IOmon = SpikeMonitor(IO)
PF_PKJmon = StateMonitor(PF_PKJ,'w', record=PF_PKJ['i!=j'])
@network_operation(dt=10*ms)
def update_active():
    GRindex = GRmon.i[GRmon.t > defaultclock.t-window]  # index for the active neuron
    IOindex = IOmon.i[IOmon.t > defaultclock.t-window]
    bothIndex = intersect1d(array(GRindex),array(IOindex))
    synapseIndexArray = PF_PKJ.i[0:np.size(PF_PKJ.i)]
    synapseIndexMat = []
    for j in range(0,np.size(bothIndex)):
        if np.size(bothIndex) > 0:
            synapseIndex = np.where(synapseIndexArray==bothIndex[j])[0]
            synapseIndexMat.extend(synapseIndex)    
    print synapseIndexMat
    PF_PKJ.w_[synapseIndexMat] *= 0.995 
mon = SpikeMonitor(PKJ)
trace = StateMonitor(PKJ, 'v', record=True)
trace2 = StateMonitor(IO,'v', record=True)
trace3 = StateMonitor(GR, 'v', record=True)

GR.I = 0*nA
IO.I[5] = 5*nA
IO.I[6] = 5*nA
IO.I[9] = 5*nA
GR.I = 5*nA
#IO.I = 100*nA
run(1000*ms,report='text')

figure(1)
plot(trace.t/ms,trace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('PKJ Voltage')

figure(2)
plot(trace2.t/ms,trace2.v[0]/mV)
title('IO Voltage')

figure(3)
plot(trace3.t/ms,trace3.v[0]/mV)

