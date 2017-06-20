# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:40:05 2016

@author: Elijah
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 16:39:54 2016

@author: Elijah
"""
# Iz neuron network

from brian2 import *
start_scope()

#%%
####################################################################
############ variables #############################################
####################################################################


n = 100;
duration = 200*ms
taus = 1*ms
neuronType1 ='tc'
neuronType2 = 'msn'

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

#%%
####################################################################
############ morphology ############################################
####################################################################

#morpho = Soma(13*um) # what you would expect in mouse msn
#morpho.L = Cylinder(diameter=1*um,length=500*um,n=5) # these are for components besides soma... probably too much detail
#morpho.LL = Cylinder(diameter=1*um,length=20*um,n=2)
#morph.R = Cylinder(diameter=1*um,length=100*um,n=5)




#%%
####################################################################
############ eqs ###################################################
####################################################################


eqs =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a*(b*v-u) : volt/second
I : amp
'''

eqs2 =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a2*(b2*v-u) : volt/second
I : amp
'''

if neuronType1 == 'msn':
    eqs = '''
    dv/dt = (1/ms/mV)*v**2+(20/ms)*v+80*mV/ms-u+I/nF : volt
    du/dt = a*(b*v-u) : volt/second
    I : amp
    '''
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2791037/

reset='''
v=c
u += d
'''
reset2='''
v=c2
u +=d2
'''
#%%
###############################################################################
############ neurongroups #####################################################
###############################################################################


# group 
group=NeuronGroup(n,eqs,threshold='v > 30*mV',
                    reset=reset,
                    method='euler')
group.v = c
group.u = b*c

# group2 
group2=NeuronGroup(n,eqs2,threshold='v > 30*mV',
                    reset=reset2,
                    method='euler')
group2.v = c2
group2.u = b2*c2

# poisson group
P = PoissonGroup(n, np.arange(100)*Hz + 10*Hz)

#%%
####################################################################
############ synapses ##############################################
####################################################################


S = Synapses(group,group2,on_pre='v+=1*mV')
S.connect(condition='i!=j',p=0.5)
S.delay = 2*ms

S2 = Synapses(P, group, on_pre='v+=1*mV')
S2.connect(condition='i!=j',p=0.2)

S3 = Synapses(group2,group,on_pre='v+=-150*mV')
S3.connect(condition='i!=j',p=0.2)
S3.delay = 5*ms

#%%
####################################################################
############ monitors ##############################################
####################################################################


Spike = SpikeMonitor(group)
trace = StateMonitor(group, 'v', record=True)
Spike2 = SpikeMonitor(group2)
trace2 = StateMonitor(group2, ('v'), record=True)


#%%
####################################################################
############ run ###################################################
####################################################################
group.I = 0*nA
run(duration,report='text')
group.I = 10*nA
run(duration, report='text')
group.I = 0*nA
run(duration, report='text')

#%%
####################################################################
############ plotting ##############################################
####################################################################
figure(1)
plot(trace.v[0])
xlabel('Time')
ylabel('Voltage(v)')

figure(2)
plot(trace2.v[0])
xlabel('Time')
ylabel('Voltage(v)')