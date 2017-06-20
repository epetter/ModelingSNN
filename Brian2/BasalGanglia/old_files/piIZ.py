# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 16:39:54 2016

@author: Elijah
"""
# Iz neurons connected

from brian2 import *

start_scope()

n = 100;
duration = 200*ms
#taus = 1*ms
neuronType1 ='ch'
neuronType2 = 'lts'

if neuronType1 == 'rs':
    c = -66*mV
    b= 0.2/ms
    a = 0.02/ms
    d = 8*mV/ms
    
if neuronType1 == 'ib':
    c = -55*mV
    b= 0.2/ms
    a = 0.02/ms
    d = 4*mV/ms

if neuronType1 == 'ch':
    c = -50*mV
    b= 0.2/ms
    a = 0.02/ms
    d = 2*mV/ms
    
if neuronType1 == 'fs':
    c = -66*mV
    b= 0.2/ms
    a = 0.1/ms
    d = 2*mV/ms    
    
if neuronType1 == 'tc':
   c = -63*mV
   b = 0.25/ms
   a = 0.02/ms
   d = 0.05*mV/ms
   
if neuronType1 == 'rz':
   c = -66*mV
   b = 0.25/ms
   a = 0.11/ms
   d = 2*mV/ms

if neuronType1 == 'lts':
   c = -66*mV
   b = 0.25/ms
   a = 0.02/ms
   d = 2*mV/ms

if neuronType1 == 'msn':   
   c = -55*mV
   b = -.20/ms
   a = 0.01/ms
   d = 15*mV/ms  
   #https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=128818#tabs-2

if neuronType2 == 'rs':
    c2 = -66*mV
    b2= 0.2/ms
    a2 = 0.02/ms
    d2 = 8*mV/ms
    
if neuronType2 == 'ib':
    c2 = -55*mV
    b2 = 0.2/ms
    a2 = 0.02/ms
    d2 = 4*mV/ms

if neuronType2 == 'ch':
    c2 = -50*mV
    b2= 0.2/ms
    a2 = 0.02/ms
    d2 = 2*mV/ms
    
if neuronType2 == 'fs':
    c2 = -66*mV
    b2= 0.2/ms
    a2 = 0.1/ms
    d2 = 2*mV/ms    
    
if neuronType2 == 'tc':
   c2 = -63*mV
   b2 = 0.25/ms
   a2 = 0.02/ms
   d2 = 0.05*mV/ms
   
if neuronType2 == 'rz':
   c2 = -66*mV
   b2 = 0.25/ms
   a2 = 0.11/ms
   d2 = 2*mV/ms

if neuronType2 == 'lts':
   c2 = -66*mV
   b2 = 0.25/ms
   a2 = 0.02/ms
   d2 = 2*mV/ms

if neuronType2 == 'msn':   
   c2 = -55*mV
   b2 = -.20/ms
   a2 = 0.01/ms
   d2 = 15*mV/ms  
   #https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=128818#tabs-2

eqs =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a*(b*v-u) : volt/second
I : amp
'''
reset='''
v=c
u += d
'''
#eqs = 
P = PoissonGroup(n, np.arange(100)*Hz + 10*Hz)
G = NeuronGroup(n, eqs,threshold='v > 30*mV',
                    reset=reset,
                    method='euler')
G.v = c
G.u = b*c


S = Synapses(P, G, on_pre='v+=10*mV')
S.connect(condition='i!=j',p=0.2)
Sp = StateMonitor(G, 'v', record=True)
run(20*ms)

figure()
plot(Sp.v[0])

neverWork = 0
if neverWork == 1:

    #eqs = 
    P = PoissonGroup(100, np.arange(100)*Hz + 10*Hz)
    G = NeuronGroup(100, 'dv/dt = -v / (10*ms) : volt')
    #G.v = c
    #G.u = b*c
    
    
    S = Synapses(P, G, on_pre='v+=0.1*mV')
    S.connect(j='i')
    Sp = StateMonitor(G, 'v', record=True)
    run(20*ms)
    
    figure()
    plot(Sp.v[0])
