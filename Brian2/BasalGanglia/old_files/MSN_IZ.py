# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 13:31:20 2016

@author: elijah
"""

from brian2 import *

# IZ MSN
start_scope()
duration = 100*ms

#%% Equations

eqs = '''
dv/dt = (((v+80)*(v+25)-u+I)/ms)/50 : 1
du/dt = 0.01*(-20*(v+80)-u)/ms : 1
I : 1
'''

reset='''
v=-55
u += 150
'''

#%% Neurons
P = NeuronGroup(10, model=eqs,
threshold='v>40', reset=reset, method='euler')

P.v = -55

#%% Neurons


#%% Monitor
S = SpikeMonitor(P)
Trace = StateMonitor(P,'v',record=True)

#%% Run
P.I = 400
run(duration)
P.I = 0
run(duration)
P.I = 430
run(duration)


#%% Plot

figure(1)
plot(Trace.t/ms,Trace.v[0])
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('MSN Voltage')