# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 18:47:28 2016

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
neuronType2 = 'rs'

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

c3 = -66*mV
b3= 0.2/ms
a3 = 0.1/ms
d3 = 2*mV/ms


#%%
####################################################################
############ eqs ###################################################
####################################################################


eqs =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a*(b*v-u) : volt/second
I : amp
x : metre
'''

eqs2 =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a2*(b2*v-u) : volt/second
I : amp
x : metre
'''

eqs3 = '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a3*(b3*v-u) : volt/second
I : amp
x : metre
'''

if neuronType1 == 'msn':
    eqs = '''
    dv/dt = (1/ms/mV)*v**2+(20/ms)*v+80*mV/ms-u+I/nF : volt
    du/dt = a*(b*v-u) : volt/second
    I : amp
    '''
    
if neuronType2 == 'msn':
    eqs2 = '''
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
reset3='''
v=c3
u +=d3
'''

neuron_spacing1 = 50*umetre
width = n/4.0*neuron_spacing1

neuron_spacing2 = 500*umetre
width = n/4.0*neuron_spacing2

neuron_spacing3 = 20*umetre
width = n/4.0*neuron_spacing3

neuron_spacing4 = 20*umetre
width = n/4.0*neuron_spacing4
#%%
###############################################################################
############ neurongroups #####################################################
###############################################################################


# thalamo-cortical cells
group=NeuronGroup(n, eqs,threshold='v > 30*mV',
                    reset=reset, 
                    method='euler')
group.v = c
group.u = b*c
group.x = 'i*neuron_spacing1'

# medium spiny neurons 
group2=NeuronGroup(n,eqs2,threshold='v > 40*mV',
                    reset=reset2,
                    method='euler')
group2.v = c2
group2.u = b2*c2
group.x = 'i*neuron_spacing2'

# SNr neurons
group3=NeuronGroup(n,eqs3,threshold='v>30*mV', reset=reset3,method='euler')

group3.v = c3
group3.u = b3*c3
group.x = 'i*neuron_spacing3'

# thalamic neurons
group4 =NeuronGroup(n,eqs,threshold='v>30*mV',reset=reset,method='euler')

group4.v = c
group4.u = b*c
group.x = 'i*neuron_spacing4'

# poisson groups
P = PoissonGroup(n, np.arange(100)*Hz + 10*Hz) # input to cortical neurons
P2 = PoissonGroup(n,np.arange(100)*Hz + 10*Hz) # input to SNr neurons
P3 = PoissonGroup(n,np.arange(100)*Hz + 10*Hz) # input to thalamic neurons

# SNr/ GPi neurons
#group3 = NeuronGroup(n,eqs, threshold='v> 30*mV)

#%%
####################################################################
############ synapses ##############################################
####################################################################


#S = Synapses(group,group2,'w : 1', on_pre='v += w*mV')
#S.connect(condition='i!=j',p=0.5)
#S.delay = 2*ms
#S.w = '(exp(-(x_pre-x_post)**2/(2*width**2)))' # input varies with distance
taupre = taupost = 20*ms
wmax = 0.01
Apre = 0.01
Apost = -Apre*taupre/taupost*1.05

S = Synapses(group,group2,             
             '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             ''',
             on_pre='''
             v_post += w*mV
             apre += Apre
             w = clip(w+apost, 0, wmax)
             ''',
             on_post='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             ''')
S.connect(condition='i!=j',p=0.5)
S.delay = 2*ms
S.w = '(exp(-(x_pre-x_post)**2/(2*width**2)))' # input varies with distance



S2 = Synapses(P, group, on_pre='v+=1*mV')
S2.connect(condition='i!=j',p=0.2)
S2.delay = 2*ms

S3 = Synapses(group2,group3,             
             '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             ''',
             on_pre='''
             v_post += w*mV
             apre += Apre
             w = clip(w+apost, 0, wmax)
             ''',
             on_post='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             ''')
S3.connect(condition='i!=j',p=0.2)
S3.delay = 5*ms
S3.w = '(10*exp(-(x_pre-x_post)**2/(2*width**2)))' # input varies with distance

S4 = Synapses(group3,group4,
              '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             ''',
             on_pre='''
             v_post += w*mV
             apre += Apre
             w = clip(w+apost, 0, wmax)
             ''',
             on_post='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             ''')
S4.connect(condition='i!=j',p=0.3)
S4.delay = 4*ms
S4.w = '(exp(-(x_pre-x_post)**2/(2*width**2)))' # input varies with distance

S5 = Synapses(P2,group3,on_pre='v+=.75*mV')
S5.connect(condition='i!=j',p=1)
S5.delay = 2*ms

S6 = Synapses(group4,group,
              '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             ''',
             on_pre='''
             v_post += w*mV
             apre += Apre
             w = clip(w+apost, 0, wmax)
             ''',
             on_post='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             ''')
S6.connect(condition='i!=j',p=0.2)
S6.delay = 2*ms

S7 = Synapses(P3,group4,on_pre='v+=.75*mV')
S7.connect(condition='i!=j',p=1)
S7.delay = 2*ms

# i_x i_y can give you some distance... delay can be a function of distance
# lfp np.sum(Im/r). Sum total current divided by distance from a location

#%%
####################################################################
############ monitors ##############################################
####################################################################


Spike = SpikeMonitor(group)
trace = StateMonitor(group, 'v', record=True)
Spike2 = SpikeMonitor(group2)
trace2 = StateMonitor(group2, ('v'), record=True)
Spike3 = SpikeMonitor(group3)
trace3 = StateMonitor(group3, ('v'), record=True)
Spike4 = SpikeMonitor(group4)
trace4 = StateMonitor(group4, ('v'), record=True)

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
plot(trace.t/ms,trace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')

figure(2)
plot(trace2.t/ms,trace2.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')

figure(3)
plot(trace3.t/ms,trace3.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')

figure(4)
plot(trace4.t/ms,trace4.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')