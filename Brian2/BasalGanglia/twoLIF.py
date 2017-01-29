# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 16:37:02 2016

@author: Elijah
"""

from brian2 import *
import random
import numpy as np 
import matplotlib.pyplot as plt
import time
t0 = time.clock()

#this is a dimensional form of Izhi neurons

defaultclock.dt =0.01*ms
div=defaultclock.dt
duration = 1500*ms
area=20000*umetre**2
El=-10*mV
gl=0.1*msiemens/cm**2
Cm=1.0*ufarad/cm**2
tau_ampa=2*ms


eqs_il = '''
il = gl * (El-v) :amp/meter**2
'''

eqs = '''
dv/dt = (il +I/area + (g_ampa*(60*mV-v))/area)/Cm:  volt
dg_ampa/dt=-g_ampa/tau_ampa : siemens
I: amp
'''
eqs += (eqs_il) 




Grs1 = NeuronGroup(1,  eqs,
                    threshold='v > 40*mV',
                    reset='v = El',
                    clock=Clock(defaultclock.dt),method='euler')


Grs2 = NeuronGroup(1,  eqs,
                    threshold='v > 40*mV',
                    reset='v = El',
                    clock=Clock(defaultclock.dt),method='euler')
                    
Grs1.v=El
Grs2.v=El                    


Sr = Synapses(Grs1, Grs2, clock=Grs1.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa += w
		''')
                            
                                                                                                   
Sr.connect(i=[0],j=[0])
Sr.w=240*nS
Sr.delay=5*ms #introduces a fixed delay between the firing of the pre cell  and the postsynaptic response

M = StateMonitor(Grs1,('v','g_ampa'),record=True)
M2 = StateMonitor(Grs2,('v','g_ampa'),record=True)

Grs1.I=0*nA
Grs2.I=0*nA  
run(10*ms,report='text')
Grs1.I=20*nA
Grs2.I=0*nA 
run(1*ms,report='text')
Grs1.I=0*nA
Grs2.I=0*nA  
run(30*ms,report='text')


subplot(4,1,1)
plot(M.t/ms, M.v[0]/mV)
xlim(0, 50)
subplot(4,1,2)
plot(M.t/ms, M.g_ampa[0]/nsiemens,'r')
xlim(0, 50)
subplot(4,1,3)
plot(M2.t/ms, M2.v[0]/mV)
xlim(0, 50)
subplot(4,1,4)
plot(M2.t/ms, M2.g_ampa[0]/nsiemens,)
xlim(0, 50)
print time.clock() - t0, "seconds process time"
show()