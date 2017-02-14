# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 09:42:16 2016

@author: Elijah
"""

from brian2 import *
duration = 10000*ms

#%%
# Parameters
# thalamo-cortical Iz
c = -66*mV
b= 0.2/ms
a = 0.02/ms
d = 8*mV/ms

#%%
#Eqs
eqs =  '''
dv/dt = (((0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF)) : volt #+x)*((1./taum)/ms) : volt
du/dt = a*(b*v-u) : volt/second
dx/dt = (-x+y)*(1./taupsp) : volt
dy/dt = -y*(1./taupsp)+25.27*mV/ms+(39.24*mV/ms**0.5)*xi : volt
I : amp
'''

reset='''
v=c
u += d
'''

# Neuron model parameters
#Vr = -70*mV
#Vt = -55*mV
taum = 10*ms
taupsp = 0.325*ms
weight = 4.86*mV
# Neuron model
#eqs = Equations('''
#dV/dt = (-(V-Vr)+x)*(1./taum) : volt
#dx/dt = (-x+y)*(1./taupsp) : volt
#dy/dt = -y*(1./taupsp)+25.27*mV/ms+
#(39.24*mV/ms**0.5)*xi : volt
#''')
# Neuron groups
n_groups = 10
group_size = 100
P = NeuronGroup(N=n_groups*group_size, model=eqs,
threshold='v>30*mV', reset=reset, refractory=1*ms,
method='euler')
Pinput = SpikeGeneratorGroup(85, np.arange(85),
np.random.randn(85)*1*ms + 50*ms)
# The network structure
S = Synapses(P, P, on_pre='y+=weight')
S.connect(j='k for k in range((int(i/group_size)+1)*group_size, (int(i/group_size)+2)*group_size) '
'if i<N_pre-group_size')
Sinput = Synapses(Pinput, P[:group_size], on_pre='y+=weight')
Sinput.connect()
# Record the spikes
Mgp = SpikeMonitor(P)
Minput = SpikeMonitor(Pinput)
# Setup the network, and run it
P.v = 'c + rand() * (-55*mV - c)'
run(duration)
plot(Mgp.t/ms, 1.0*Mgp.i/group_size, '.')
plot([0, duration/ms], np.arange(n_groups).repeat(2).reshape(-1, 2).T, 'k-')
ylabel('group number')
yticks(np.arange(n_groups))
xlabel('time (ms)')
show()
