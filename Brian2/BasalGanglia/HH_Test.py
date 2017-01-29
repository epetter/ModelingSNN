# -*- coding: utf-8 -*-
"""
Created on Wed Sep 07 15:27:29 2016

@author: Elijah
"""

from brian2 import *

start_scope()

num_neurons = 100
duration = .2*second

# Parameters
area = 20000*umetre**2
Cm = 1*ufarad*cm**-2
El = 10.613*mV
EK = -(12)*mV
ENa = 115*mV

gl = 0.3*msiemens/cm**2
gNa = 36*msiemens/cm**2
gK = 120*msiemens/cm**2


#The model
eqs_ina = '''
ina=gNa * m**3 * h * (ENa-v) :  amp/meter**2
dm/dt = alpham * (1-m) - betam * m : 1
dh/dt = alphah * (1-h) - betah * h : 1
alpham = (0.1/mV) * (-(v)+25*mV) / (exp((-(v)+25*mV) / (10*mV)) -1)/ms : Hz
betam = 4 * exp(-(v)/(18*mV))/ms : Hz
alphah = (0.07) * (exp((-(v)/(20*mV))))/ms : Hz
betah = 1/(exp((-(v)+(30*mV))/(10*mV)) +1)/ms :Hz 
'''

eqs_ik = '''
ik=gK * n**4 * (EK-v) : amp/meter**2
dn/dt = alphan * (1-n) - betan * n : 1
alphan = (0.01/mV) * (-(v)+10*mV) / (exp((-(v)+10*mV) / (10*mV)) - 1)/ms : Hz
betan = 0.125*exp(-(v)/(80*mV))/ms : Hz
'''

eqs_il = '''
il = gl * (El-v) : amp/meter**2
'''

eqs = '''
dv/dt = (ina+ik+il +I/area)/Cm:  volt
I : amp
'''
eqs += (eqs_ina+eqs_ik+eqs_il) 

# Threshold and refractoriness are only used for spike counting
group = NeuronGroup(num_neurons, eqs,
                    threshold='v > 40*mV',
                    refractory='v > 40*mV',
                    method='exponential_euler')
group.v = 0
group.m=0.0529
group.n=0.3177
group.h=0.596


spikeMon = SpikeMonitor(group)
monitor2=StateMonitor(group,'v',record=True)
group.I = 0*nA
run(5.0*ms,report='text')
group.I = 12.4*nA
run(1.0*ms, report='text')
group.I = 0*nA
run(10.0*ms)
group.I = 0*nA
run(5.0*ms,report='text')
group.I = 12.4*nA
run(1.0*ms, report='text')
group.I = 0*nA
run(10.0*ms)
#group.I = '250*nA*i/num_neurons'

run(duration)


figure(1)
plot(monitor2.t/ms, monitor2.v[0]/mV) #plot the voltage for neuron 0 (index starts at 0)
ylim(-20,120) #set axes limits
xlabel('Time')
ylabel('Voltage')

figure(2)
plot(group.I/nA, spikeMon.count/duration)
xlabel('I (nA)')
ylabel('Firing rate (sp/s)')