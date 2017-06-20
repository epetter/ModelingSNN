# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 19:38:28 2016

@author: Elijah
"""

from brian2 import *

duration = 100*ms
G=0*arange(0,duration/(0.01*ms))
G[3500]=1.0 #time in timesteps that an event takes place 3500 is 3500*.01*ms or 35.0 ms
G[7000]=1.0
A1=TimedArray(G,dt=0.01*ms)
I=A1
vt = 1.*mV         # Spiking threshold
memc = 200.0*pfarad  # Membrane capacitance
bgcurrent = 200*pA   # External current
tau_m=2*ms
tau_ampa=2*ms
g_synpk=0.3*nsiemens
transwdth=1.0*ms


eqs_neurons='''
dv/dt=-v/tau_m + I(t)*mag : volt
dg_ampa/dt = -g_ampa/tau_ampa : siemens

Tr_pre=.25*(tanh((t/ms-tspike/ms)/.005)-tanh((t/ms-(tspike/ms +transwdth/ms))/.005)):1

mag: volt/second
tspike:second
'''

# ###########################################
# Initialize neuron group
# ###########################################

neurons = NeuronGroup(1, model=eqs_neurons, threshold='v > vt',
                      reset='v=0*mV;g_ampa=g_synpk;tspike=t',refractory='0.3*ms')

neurons.mag=150.0*mV/ms
neurons.v=0.0*mV
neurons.tspike=-100.0*ms  #needed so a spike does not happen at time 0
# Comment these two lines out to see what happens without Synapses


M = StateMonitor(neurons, ('v','g_ampa','Tr_pre'), record=True)
sm = SpikeMonitor(neurons)

run(duration)
figure(1)
subplot(3,1,1)
plot(M.t/ms, M.v[0]/mV, '-b')
xlim(0,duration/ms)
subplot(3,1,2)
plot(M.t/ms, M.g_ampa[0]/nsiemens, '-g', lw=2)

subplot(3,1,3)
plot(M.t/ms, M.Tr_pre[0], '-r', lw=2)
xlim(0,duration/ms)
show()