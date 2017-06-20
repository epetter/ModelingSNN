# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 09:42:16 2016

@author: Elijah
"""

from brian2 import *
duration = 1000*ms

#%%
# Parameters
# thalamo-cortical Iz
c = -66#-66*mV
b= 0.2#0.2/ms
a = 0.02#0.02/ms
d = 8#8*mV/ms

#%%

taum = 10*ms
taupsp = 0.325*ms
weight = 4.86*mV


# Missing the (1/.taum)
#eqs =  '''
#dv/dt = (((0.04)*v**2+(5)*v+140-u+I)+x)/ms : 1
#du/dt = a*(b*v-u)/ms : 1
#dx/dt = ((-x+y/mV)*(1./taupsp)) : 1
#dy/dt = -y*(1./taupsp)+25.27*mV/ms+(39.24*mV/ms**0.5)*xi : volt
#I : 1
#'''

# equations similar to original brian synfire example... 
# The groups are limited to about 4 and span larger distances 
eqs =  '''
dv/dt = ((((0.04)*v**2+(5)*v+140-u+I)+x)*((1./taum)*ms))/ms : 1
du/dt = a*(b*v-u)/ms : 1
dx/dt = ((-x+y/mV)*(1./taupsp)) : 1
dy/dt = -y*(1./taupsp)+25.27*mV/ms+(39.24*mV/ms**0.5)*xi : volt
I : 1
'''

MSNeqs = '''
dv/dt = (((v+80)*(v+25)-u+I)/ms)/50 : 1
du/dt = 0.01*(-20*(v+80)-u)/ms : 1
I : 1
'''

reset='''
v=c
u += d
'''

MSNreset='''
v=-55
u += 150
'''

# Neuron model parameters
#Vr = -70*mV
#Vt = -55*mV

# Neuron model
#eqs = Equations('''
#dV/dt = (-(V-Vr)+x)*(1./taum) : volt
#dx/dt = (-x+y)*(1./taupsp) : volt
#dy/dt = -y*(1./taupsp)+25.27*mV/ms+
#(39.24*mV/ms**0.5)*xi : volt
#'''
#%% Neuron groups
n_groups = 10
group_size = 100

P = NeuronGroup(N=n_groups*group_size, model=eqs,
threshold='v>30', reset=reset, refractory=1*ms,
method='euler')

Pinput = SpikeGeneratorGroup(85, np.arange(85),
np.random.randn(85)*1*ms + 50*ms)

MSN = NeuronGroup(2,model=MSNeqs,threshold='v>40',reset=reset,refractory=2*ms,method='euler')


#%% The network structure
S = Synapses(P, P, on_pre='y+=weight')
S.connect(j='k for k in range((int(i/group_size)+1)*group_size, (int(i/group_size)+2)*group_size) '
'if i<N_pre-group_size')
Sinput = Synapses(Pinput, P[:group_size], on_pre='y+=weight')
Sinput.connect()

Cortex_MSN = Synapses(P,MSN,on_pre='v+=6')
Cortex_MSN.connect(condition='i!=j',p=0.1)

#%% Record the spikes
Mgp = SpikeMonitor(P)
#Minput = SpikeMonitor(Pinput)
# Setup the network, and run it
P.v = 'c + rand() * (-55 - c)'

MSNspikes = SpikeMonitor(MSN)
MSNtrace = StateMonitor(MSN,'v',record=True)

#%%
# run 

run(duration)

#%% Plotz
plot(Mgp.t/ms, 1.0*Mgp.i/group_size, '.')
plot([0, duration/ms], np.arange(n_groups).repeat(2).reshape(-1, 2).T, 'k-')
ylabel('group number')
yticks(np.arange(n_groups))
xlabel('time (ms)')
show()

plot(MSNtrace.t/ms,MSNtrace.v[0])
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('MSN Voltage')