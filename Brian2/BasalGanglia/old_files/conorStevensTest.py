# -*- coding: utf-8 -*-
"""
Created on Wed Sep 07 18:06:10 2016

@author: Elijah
"""

from brian2 import *

start_scope()

num_neurons = 100
duration = 2*second

# Parameters
area = 20000*umetre**2
Cm = 1*ufarad*cm**-2
El = -17.0*mV
EK = -72*mV
ENa = 55.0*mV
EA= -75.0*mV

#gl = 0.3*msiemens/cm**2
#gNa = 120*msiemens/cm**2
#gK = 20.0*msiemens/cm**2
#gA=47.7*msiemens/cm**2

gl = 0.3*msiemens/cm**2
gNa = 120.0*msiemens/cm**2
gK = 20*msiemens/cm**2
gA=  47.7*msiemens/cm**2

#The model
eqs_ina = '''
ina=gNa * m**3 * h * (ENa-v) :  amp/meter**2
dm/dt = alpham * (1-m) - betam * m : 1
dh/dt = alphah * (1-h) - betah * h : 1
alpham = 0.38/mV*(v+29.7*mV)/(1-exp(-0.1*(v+29.7*mV)/mV ) )/ms : Hz
betam = 15.2*exp(-0.0556*(v+54.7*mV)/mV)/ms : Hz
alphah = 0.266*exp(-0.05*(v+48*mV)/mV)/ms : Hz
betah = 3.8/(1+exp(-0.1*(v+18.*mV)/mV))/ms : Hz
'''


eqs_iA = '''
iA=gA * a**3 * b * (EA-v) :  amp/meter**2
da/dt = alphaa * (1-a) - betaa * a : 1
db/dt = alphab * (1-b) - betab * b : 1
alphaa = (0.0761*exp(0.0314*(v+94.22*mV)/mV)/(1+exp(0.0346*(v+1.17*mV)/mV))) ** (1/3)/ms : Hz 
betaa = ((0.3632*mV)/mV + 1.158/(1 + exp(0.0497*(v + 55.96*mV)/mV)))**-1/ms : Hz
alphab = (1/(1+exp(0.0688*(v+53.3*mV)/mV))) ** 4/ms : Hz
betab = ((1.24*mV)/mV + 2.678/(1 + exp(0.0624*(v + 50*mV)/mV)))**-1/ms : Hz
'''

#iA=gA * a**3 * b * (EA-v) :  amp/meter**2
#da/dt = alphaa * (1-a) - betaa * a : 1
#db/dt = alphab * (1-b) - betab * b : 1
#alphaa = ((0.0761*exp(0.0314*(v+94.22*mV)/mV))/(1+exp(0.0346*(v+1.17*mV)/mV)) ** (1/3))/ms : Hz 
#betaa = ((0.3632*mV)/mV + 1.158/(1 + exp(0.0497*(v + 55.96*mV)/mV)))/ms : Hz
#alphab = (1/(1+exp(0.0688*(v+55.3*mV)/mV)) ** 4)/ms : Hz
#betab = ((1.24*mV)/mV + 2.678/(1 + exp(0.0624*(v + 50*mV)/mV)))/ms : Hz
#'''

eqs_ik = '''
ik=gK * n**4 * (EK-v):amp/meter**2
dn/dt = alphan * (1-n) - betan * n : 1
alphan = (0.02*(v+45.7*mV)/mV)/(1-exp(-0.1*(v+45.7*mV)/mV))/ms : Hz
betan = 0.25*exp(-0.0125*(v+55.7*mV)/mV)/ms : Hz

'''
eqs_il = '''
il = gl * (El-v) :amp/meter**2
'''

eqs = '''
dv/dt = (ina+ik+il+iA+I/area)/Cm:  volt
I : amp
'''
eqs += (eqs_ina+eqs_ik+eqs_il+eqs_iA) 

# Threshold and refractoriness are only used for spike counting
group = NeuronGroup(num_neurons, eqs,
                    threshold='v > 40*mV',
                    refractory='v > 40*mV',
                    method='exponential_euler')
group.v = -160.0*mV
group.m=0.0529
group.n=0.3177
group.h=0.596

monitor3=StateMonitor(group,'iA',record=True)
monitor4=StateMonitor(group, 'ina',record=True)

monitor2=StateMonitor(group,'v',record=True)
s = SpikeMonitor(group)
group.I = 0*nA
run(25.0*ms,report='text')
group.I = 1.7*nA
run(25.0*ms, report='text')
group.I = 0*nA
run(25.0*ms)
#group.I = 1.7*nA
#3run(25*ms)


#group.I = '(7*nA * i) / num_neurons'

#monitor = SpikeMonitor(group)

#run(duration)

figure(1)
plot(group.I/nA, monitor.count/duration)
xlabel('I (nA)')
ylabel('Firing rate (sp/s)')

figure(2)
plot(monitor2.t/ms, monitor2.v[0]/mV) #plot the voltage for neuron 0 (index starts at 0)
ylim(-80,60) #set axes limits
xlabel('Time(ms)')
ylabel('voltage(mV)')

#figure(3)
#plot(monitor3.t/ms, monitor3.iA[0]/(amp/meter**2))

#figure(4)
#plot(monitor2.t/ms, monitor2.v[0]/mV)
#plot(monitor3.t/ms, monitor3.iA[0]/(amp/meter**2))

#figure(5) 
#plot(monitor4.t/ms, monitor4.ina[0]/(amp/meter**2))
#show()

#figure(6)
#plot(group.I/nA,monitor.count/duration)