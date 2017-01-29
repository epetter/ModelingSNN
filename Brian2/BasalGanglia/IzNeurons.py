# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 14:21:33 2016

@author: Elijah
"""

# Iz neurons

from brian2 import *

start_scope()

n = 1;
duration = 200*ms
neuronType ='msn'

if neuronType == 'rs':
    c = -66*mV
    b= 0.2/ms
    a = 0.02/ms
    d = 8*mV/ms
    
if neuronType == 'ib':
    c = -55*mV
    b= 0.2/ms
    a = 0.02/ms
    d = 4*mV/ms

if neuronType == 'ch':
    c = -50*mV
    b= 0.2/ms
    a = 0.02/ms
    d = 2*mV/ms
    
if neuronType == 'fs':
    c = -66*mV
    b= 0.2/ms
    a = 0.1/ms
    d = 2*mV/ms    
    
if neuronType == 'tc':
   c = -63*mV
   b = 0.25/ms
   a = 0.02/ms
   d = 0.05*mV/ms
   
if neuronType == 'rz':
   c = -66*mV
   b = 0.25/ms
   a = 0.11/ms
   d = 2*mV/ms

if neuronType == 'lts':
   c = -66*mV
   b = 0.25/ms
   a = 0.02/ms
   d = 2*mV/ms

if neuronType == 'msn':   
   c = -55*mV
   b = -.20/ms
   a = 0.01/ms
   d = 15*mV/ms  
   #https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=128818#tabs-2



eqs =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a*(b*v-u) : volt/second
I : amp
'''

if neuronType == 'msn':
   eqs = ''' # Izikevich 2007 book 
   dv/dt = (((v+80)*(v+25)-u+I)/ms)/50 : 1
   du/dt = 0.01*(-20*(v+150))/ms : 1
   I : 1
'''
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2791037/

reset='''
v=c
u += d
'''


group=NeuronGroup(n,eqs,threshold='v > 30*mV',
                    reset=reset,
                    method='euler')
group.v = c
group.u = b*c



S = SpikeMonitor(group)
trace = StateMonitor(group, 'v', record=True)
group.I = 0*nA
run(duration,report='text')
group.I = 10*nA
run(duration, report='text')
#group.I = 0*nA
#run(duration, report='text')
#group.I = .75*nA
#run(duration, report='text')
#group.I = '250*nA*i/n'
#run(duration,report='text')

figure(1)
plot(trace.v[0])
xlabel('Time')
ylabel('Voltage(v)')