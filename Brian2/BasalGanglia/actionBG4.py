# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:21:53 2017

@author: elijah
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 11:03:05 2016

@author: elijah
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 10:05:50 2016

@author: elijah
"""

###############
# Basal Ganglia for action selection
# EP 170316
###############

# Add stochasticity xi*second**.5

# make cortical activity more random across channels

#%% import and setup 
from brian2 import *
import numpy as np
import scipy as sp
import scipy.stats
import matplotlib.pyplot as plt
#prefs.codegen.target = 'weave'
start_scope()

# options
synfire = 0
SNr_collaterals = 0
recordz = 0
plotz = 0
action_thresh = 10

# Tests/ experiments to run 
sequence = 0
popFiring = 0
cortex_D1_action = 0 # a test to see if increased Cortex D1 strength can choose an action 
learnAction = 1 # test to see if an action can be learned  
test_DA = 0 # test to look at DA firing 

# variables 
pop_duration = 11000*ms # the duration to run simulations for population firing rates. This was 11 seconds in Humphries et al., 2006; 
sequence_duration = 1500*ms # As there are three stages this will result in a 3 seconds simulation
learn_duration = 500000*ms 
synfire_duration = 100*ms # a quick test to make sure the synfire chain is functioning correctly 
cortex_D1_duration = 3000*ms # a test of whether or not I can achieve more actions just through cortical-D1 plasticity 
DA_duration = 100*ms
binSize = 100*ms 

# cortical
taum = 5*ms
taupsp = 0.325*ms#0.325*ms
weight = 4.86*mV

# Cortical-MSN plasticity
MSN_High = -60
SNr_thresh = 25*Hz 

#%%
# Parameters

n = 10 # number of neurons
w2 = n # for connectivity
numChannels = 3 # the number of action channels... to scale STN and GPe 

# TC Neuron
c0 = -63
b0 = 0.25
a0 = 0.02
d0 = 0.05

# RS Iz
c = -66#-66*mV
b= 0.2#0.2/ms 
a = 0.02#0.02/ms
d = 8#8*mV/ms

# FS Iz
c2 = -66
b2= 0.2
a2 = 0.1
d2 = 2  

# STDP variables
taupre = taupost = 20*ms
wmax = 50 # weights can go from 0 to 1
Apre = 0.01
Apost = -Apre*taupre/taupost*1.05
wScale = 1 # amount to scale weights by

# MSN plasticity
traceTau = traceTauPost = 20*ms  # this will change integration window # was 1200
traceConstant = 0.01 # this will change maximum weight change
tracePostConstant = -traceConstant*traceTau/traceTauPost*1.05

# Starting Weights... Different weights based upon synaptic distributions Humphries, et al., 2006 
s = 10 # Somatic 
p = 8 # proximal
d = 5 # distal 

# firing rate window
integration_window = 100*ms

#%% Equations 
############ Neuron Eqs

# equations similar to original brian synfire example... 
# The groups are limited to about 4 and span larger distances 
synfire_eqs =   '''
dv/dt = ((((0.04)*v**2+(5)*v+140-u+I)+x)*((1./taum)*ms))/ms : 1
du/dt = a*(b*v-u)/ms : 1
dx/dt = ((-x+y/mV)*(1./taupsp)) : 1
dy/dt = -y*(1./taupsp)+25.27*mV/ms+(39.24*mV/ms**0.5)*xi : volt
I : 1
'''

eqs = '''
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = a0*(b0*v-u)/ms : 1
I : 1
''' 
eqs2 = '''
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = a*(b*v-u)/ms : 1
I : 1
'''
eqs3= '''
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = a2*(b2*v-u)/ms : 1
I : 1
''' 

GPeEqs=''' # Mandali et al., 2015
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = 0.1*(0.2*v-u)/ms : 1 
I : 1
''' 
    
MSNeqs = ''' # Izikevich 2007 book 
dv/dt = (((v+80)*(v+25)-u+I)/ms)/50 : 1
du/dt = 0.01*(-20*(v+80)-u)/ms : 1
I : 1
'''

STNeqs = ''' # Mandali et al., 2015
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = 0.005*(0.265*v-u)/ms : 1 
I : 1
'''  
 
SNrEqs = ''' # Mandali et al., 2015
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = 0.1*(0.2*v-u)/ms : 1 
I : 1
'''   

############ Resets 
reset0='''
v=c0
u += d0
''' 

reset='''
v=c
u += d
'''

reset3='''
v=c2
u += d2
'''

GPeReset='''
v=(-65)
u+=2
'''

MSNreset = '''
v = -55
u+=150
'''

STNreset='''
v=-65
u+=1.5
'''

SNrReset='''
v=-65
u+=2
'''

############ Plasticity Equations 
# weight equations
weightEqs='''
                w : 1
'''

subW='''
             v -= w
             '''

addW='''
             v += w
             '''

Cortex_MSN='''
             v += w
             '''

# Cortical Striatal              
MSNstdp=  '''
             w : 1
             dtraceCon/dt = -traceCon/traceTau : 1 (event-driven)
             dtraceConPost/dt = -traceConPost/traceTauPost : 1 (event-driven)
             
             '''    
MSNpre='''
             v_post += w
             traceCon += traceConstant
             '''

MSNpost ='''
             traceConPost += tracePostConstant
             
             '''

# DA plasticity
DAstdp = '''
             w : 1
             dtraceCon/dt = -traceCon/traceTau : 1 (event-driven)
             dtraceConPost/dt = -traceConPost/traceTauPost : 1 (event-driven)
             '''   
onpreDA_D1='''
             v_post += w*wScale
             traceCon += traceConstant
             '''                                              
onpostDA='''
             traceConPost += tracePostConstant 
         '''
onpreDA_D2='''
             v_post -= w
             traceCon += traceConstant
             '''            
#%% Neuron Groups 
             
if synfire == 1:
   n_groups = 10
   group_size = 80 # This will result in 800 neurons which is what Laje and Buonomano 2013 used in their RNN
   w3 = n_groups*group_size             

############ Poisson Neurons 
CortexL_Poisson = PoissonGroup(n, np.arange(n)*Hz + 10*Hz) # input to cortical neurons
CortexNL_Poisson = PoissonGroup(n, np.arange(n)*Hz + 10*Hz) # input to cortical neurons
CortexWL_Poisson = PoissonGroup(n, np.arange(n)*Hz + 10*Hz) # input to cortical neurons

############ WTA neurons 
LeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
NoLeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
WrongLeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
#InhInter = NeuronGroup(5,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterNL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterWL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')

############ Cortical Neurons 
if synfire != 1:
   CortexL = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset0,method='euler')
    
   CortexL.v = c0
   CortexL.u = b0*c0
    
   CortexNL = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset0,method='euler')
    
   CortexNL.v = c0
   CortexNL.u = b0*c0

   CortexWL = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset0,method='euler')

   CortexWL.v = c0
   CortexWL.u = b0*c0
   
elif synfire == 1:
     CortexL = NeuronGroup(N=n_groups*group_size, model=synfire_eqs,
                           threshold='v>30', reset=reset0, refractory=2*ms,
                           method='euler')
     CortexL.v = c0
     CortexL.u = b0*c0                     
                             
     CortexNL = NeuronGroup(N=n_groups*group_size, model=synfire_eqs,
                           threshold='v>30', reset=reset0, refractory=2*ms,
                           method='euler')
     CortexNL.v = c0
     CortexNL.u = b0*c0                               

     CortexWL = NeuronGroup(N=n_groups*group_size, model=synfire_eqs,
                            threshold='v>30', reset=reset0, refractory=2*ms,
                            method='euler')

     CortexWL.v = c0
     CortexWL.u = b0*c0        

     # input to the cortex 
     Pinput = SpikeGeneratorGroup(85, np.arange(85),
                                 np.random.randn(85)*1*ms + 50*ms)                       
                             
############ Striatal Neurons  
# D1 
D1_L = NeuronGroup(n*2,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D1_L.v = c
D1_L.u = b*c

D1_NL = NeuronGroup(n*2,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D1_NL.v = c
D1_NL.u = b*c

D1_WL = NeuronGroup(n*2,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D1_WL.v = c
D1_WL.u = b*c

# D2
D2_L = NeuronGroup(n*2,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D2_L.v = c
D2_L.u = b*c 

D2_NL = NeuronGroup(n*2,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D2_NL.v = c
D2_NL.u = b*c 

D2_WL = NeuronGroup(n*2,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D2_WL.v = c
D2_WL.u = b*c 

############ Dopamine Neurons 
DA = NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')
DA.v = c
DA.u = b*c

############ GPe Neurons
GPe_L = NeuronGroup(n,GPeEqs,threshold='v>30',reset=GPeReset,method='euler') # Humphries, et al., 2006 

GPe_L.v = -65
GPe_L.u = 0.2*-65

GPe_NL = NeuronGroup(n,GPeEqs,threshold='v>30',reset=GPeReset,method='euler') # Humphries, et al., 2006 

GPe_NL.v = -65
GPe_NL.u = 0.2*-65

GPe_WL = NeuronGroup(n,GPeEqs,threshold='v>30',reset=GPeReset,method='euler') # Humphries, et al., 2006 

GPe_WL.v = -65
GPe_WL.u = 0.2*-65

############ SNr Neurons 
SNrL = NeuronGroup(n+10,SNrEqs,threshold='v>30',reset=SNrReset,method='euler') # Humphries, et al., 2006 

SNrL.v = -65
SNrL.u = 0.2*-65

SNrNL = NeuronGroup(n+10,SNrEqs,threshold='v>30',reset=SNrReset,method='euler') # Humphries, et al., 2006 

SNrNL.v = -65
SNrNL.u = 0.2*-65

SNrWL = NeuronGroup(n+10,SNrEqs,threshold='v>30',reset=SNrReset,method='euler') # Humphries, et al., 2006 

SNrWL.v = -65
SNrWL.u = 0.2*-65

############ STN Neurons 
STN = NeuronGroup(n*numChannels,STNeqs,threshold='v>30',reset=STNreset,method='euler') # Humphries, et al., 2006 

STN.v = (-65)
STN.u = (-65*1.5)

############ Thalamic Neurons 
ThalamusL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

ThalamusL.v = c0
ThalamusL.u = b0*c0

ThalamusNL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

ThalamusNL.v = c0
ThalamusNL.u = b0*c0

ThalamusWL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

ThalamusWL.v = c0
ThalamusWL.u = b0*c0

#%% Synapses

############ WTA circuit 
LeverPress_PressEx = Synapses(LeverPress,LeverPress,on_pre='v+=10') # excitatory with lower delays 
LeverPress_PressEx.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_PressEx.delay = 2*ms

LeverPress_InhNL = Synapses(LeverPress,InhInterNL,on_pre='v+=10')
LeverPress_InhNL.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_InhNL.delay = 2*ms

LeverPress_InhWL = Synapses(LeverPress,InhInterWL,on_pre='v+=10')
LeverPress_InhWL.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_InhWL.delay = 2*ms

# Wrong lever Press
WrongPress_WrongPressEx = Synapses(WrongLeverPress,WrongLeverPress,on_pre='v+=10')
WrongPress_WrongPressEx.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_WrongPressEx.delay = 2*ms

WrongPress_InhL = Synapses(WrongLeverPress,InhInterL,on_pre='v+=10')
WrongPress_InhL.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_InhL.delay = 2*ms

WrongPress_InhNL = Synapses(WrongLeverPress,InhInterNL,on_pre='v+=10')
WrongPress_InhNL.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_InhNL.delay = 2*ms

# No Lever Press 
NoPress_NoPressEx = Synapses(NoLeverPress,NoLeverPress,on_pre='v+=10')
NoPress_NoPressEx.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_NoPressEx.delay = 2*ms

NoPress_InhL = Synapses(NoLeverPress,InhInterL,on_pre='v+=10')
NoPress_InhL.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_InhL.delay = 2*ms

NoPress_InhWL = Synapses(NoLeverPress,InhInterWL,on_pre='v+=10')
NoPress_InhWL.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_InhWL.delay = 2*ms

# Inhibitory Interneuron
InhWLPress = Synapses(InhInterWL,LeverPress,on_pre='v=-20')
InhWLPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
InhWLPress.delay = 1*ms

InhNLPress = Synapses(InhInterNL,LeverPress,on_pre='v=-20')
InhNLPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
InhNLPress.delay = 1*ms

InhLWrongPress = Synapses(InhInterL,WrongLeverPress,on_pre='v=-20')
InhLWrongPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
InhLWrongPress.delay = 1*ms

InhNLWrongPress = Synapses(InhInterNL,WrongLeverPress,on_pre='v=-20')
InhNLWrongPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
InhNLWrongPress.delay = 1*ms

InhLNoPress = Synapses(InhInterL,NoLeverPress,on_pre='v=-20')
InhLNoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
InhLNoPress.delay = 1*ms

InhWLNoPress = Synapses(InhInterWL,NoLeverPress,on_pre='v=-20')
InhWLNoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
InhWLNoPress.delay = 1*ms


############ Cortical Projections 
# synfire connections 
if synfire == 1:
    CortexL_CortexL = Synapses(CortexL, CortexL, on_pre='y+=weight')
    CortexL_CortexL.connect(j='k for k in range((int(i/group_size)+1)*group_size, (int(i/group_size)+2)*group_size) '
    'if i<N_pre-group_size')
    
    CortexNL_CortexNL = Synapses(CortexNL, CortexNL, on_pre='y+=weight')
    CortexNL_CortexNL.connect(j='k for k in range((int(i/group_size)+1)*group_size, (int(i/group_size)+2)*group_size) '
    'if i<N_pre-group_size')    
    
    CortexWL_CortexWL = Synapses(CortexWL, CortexWL, on_pre='y+=weight')
    CortexWL_CortexWL.connect(j='k for k in range((int(i/group_size)+1)*group_size, (int(i/group_size)+2)*group_size) '
    'if i<N_pre-group_size')
    
    # synfire input into the cortex 
    #SinputL = Synapses(Pinput, CortexL[:group_size], on_pre='y+=weight')
    #SinputL.connect()
    
    #SinputNL = Synapses(Pinput, CortexNL[:group_size], on_pre='y+=weight')
    #SinputNL.connect()
    
    #SinputWL = Synapses(Pinput, CortexWL[:group_size], on_pre='y+=weight')
    #SinputWL.connect()    
    
# Cortex WTA
CortexL_Lever = Synapses(CortexL,LeverPress,on_pre='v+=10')
CortexL_Lever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) #i=[0,1,2],j=0)#j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
CortexL_Lever.delay = 15*ms

CortexNL_NoLever = Synapses(CortexNL,NoLeverPress,on_pre='v+=10')
CortexNL_NoLever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) #(i=[3,4,5],j=0)#j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
CortexNL_NoLever.delay = 15*ms

CortexWL_WrongLever = Synapses(CortexWL,WrongLeverPress,on_pre='v+=10')
CortexWL_WrongLever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) #i=[6,7,8],j=0)#j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
CortexWL_WrongLever.delay = 15*ms

# Cortex D1
if cortex_D1_action == 1:
    CortexL_D1L = Synapses(CortexL,D1_L,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
    CortexL_D1L.connect(j='k for k in range(i-w2, i+w2) if rand()<0.9', skip_if_invalid=True)
    CortexL_D1L.delay = 10*ms # Humphries, et al., 2006 
    CortexLD1start = 30
    CortexL_D1L.w = CortexLD1start #rand(len(Cortex_D1.i))
else:
    CortexL_D1L = Synapses(CortexL,D1_L,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
    CortexL_D1L.connect(j='k for k in range(i-w2, i+w2) if rand()<0.9', skip_if_invalid=True)
    CortexL_D1L.delay = 10*ms # Humphries, et al., 2006 
    CortexLD1start = s#*rand(len(CortexL_D1L.i))
    CortexL_D1L.w = CortexLD1start #rand(len(Cortex_D1.i))

CortexNL_D1NL = Synapses(CortexNL,D1_NL,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexNL_D1NL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.9', skip_if_invalid=True)
CortexNL_D1NL.delay = 10*ms # Humphries, et al., 2006 
CortexNL_D1NL.w = s#s*rand(len(CortexNL_D1NL.i))

CortexWL_D1WL = Synapses(CortexWL,D1_WL,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexWL_D1WL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.9', skip_if_invalid=True)
CortexWL_D1WL.delay = 10*ms # Humphries, et al., 2006 
CortexWL_D1WL.w = s#s*rand(len(CortexWL_D1WL.i))

# Cortex D2 
CortexL_D2 = Synapses(CortexL,D2_L,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexL_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.9', skip_if_invalid=True)
CortexL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexL_D2.w = s

CortexNL_D2 = Synapses(CortexNL,D2_NL,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexNL_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.9', skip_if_invalid=True)
CortexNL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexNL_D2.w = s

CortexWL_D2 = Synapses(CortexWL,D2_WL,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexWL_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.9', skip_if_invalid=True)
CortexWL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexWL_D2.w = s

# Cortex STN - Hyperdirect pathway 
CortexL_STN = Synapses(CortexL,STN,weightEqs,on_pre=addW)
CortexL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
CortexL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexL_STN.w = d #d was suggested by humphries, but not enough activation 

CortexNL_STN = Synapses(CortexNL,STN,weightEqs,on_pre=addW)
CortexNL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
CortexNL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexNL_STN.w = d

CortexWL_STN = Synapses(CortexWL,STN,weightEqs,on_pre=addW)
CortexWL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
CortexWL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexWL_STN.w = d

############ DopaminergicProjections 
# DA projections to D1
DA_D1L = Synapses(DA,D1_L,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1L.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
DA_D1L.delay = 5*ms
DA_D1L.w = s#rand(len(DA_D1L.i))

DA_D1NL = Synapses(DA,D1_NL,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1NL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
DA_D1NL.delay = 5*ms
DA_D1NL.w = s#rand(len(DA_D1NL.i))

DA_D1WL = Synapses(DA,D1_WL,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1WL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
DA_D1WL.delay = 5*ms
DA_D1WL.w = s#rand(len(DA_D1WL.i))


# DA projections to D1
DA_D2L = Synapses(DA,D2_L,DAstdp,on_pre=onpreDA_D2,on_post=onpostDA)
DA_D2L.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
DA_D2L.delay = 5*ms
DA_D2L.w = s#rand(len(DA_D1L.i))

DA_D2NL = Synapses(DA,D2_NL,DAstdp,on_pre=onpreDA_D2,on_post=onpostDA)
DA_D2NL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
DA_D2NL.delay = 5*ms
DA_D2NL.w = s#rand(len(DA_D1NL.i))

DA_D2WL = Synapses(DA,D2_WL,DAstdp,on_pre=onpreDA_D2,on_post=onpostDA)
DA_D2WL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
DA_D2WL.delay = 5*ms
DA_D2WL.w = s#rand(len(DA_D1WL.i))

############ Striatal Projections 
# D1 
D1_SNrL = Synapses(D1_L,SNrL,weightEqs,on_pre=subW)
D1_SNrL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
D1_SNrL.delay = 4*ms # Humphries, et al., 2006 
D1_SNrL.w = s # Humphries, et al., 2006  #rand(len(D1_SNrL.i))

D1_SNrNL = Synapses(D1_NL,SNrNL,weightEqs,on_pre=subW)
D1_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
D1_SNrNL.delay = 4*ms # Humphries, et al., 2006 
D1_SNrNL.w = s # Humphries, et al., 2006  #rand(len(D1_SNrNL.i))

D1_SNrWL = Synapses(D1_WL,SNrWL,weightEqs,on_pre=subW)
D1_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
D1_SNrWL.delay = 4*ms # Humphries, et al., 2006 
D1_SNrWL.w = s # Humphries, et al., 2006  #rand(len(D1_SNrWL.i))

# D2 
D2_L_GPe_L = Synapses(D2_L,GPe_L,weightEqs,on_pre=subW)
D2_L_GPe_L.connect(j='k for k in range(i-w2, i+w2) if rand()<0.8', skip_if_invalid=True)
D2_L_GPe_L.delay = 5*ms # Humphries, et al., 2006 
D2_L_GPe_L.w = np.random.choice([s,p,d],len(D2_L_GPe_L.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

D2_NL_GPe_NL = Synapses(D2_NL,GPe_NL,weightEqs,on_pre=subW)
D2_NL_GPe_NL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.8', skip_if_invalid=True)
D2_NL_GPe_NL.delay = 5*ms # Humphries, et al., 2006 
D2_NL_GPe_NL.w = np.random.choice([s,p,d],len(D2_NL_GPe_NL.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

D2_WL_GPe_WL = Synapses(D2_WL,GPe_WL,weightEqs,on_pre=subW)
D2_WL_GPe_WL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.8', skip_if_invalid=True)
D2_WL_GPe_WL.delay = 5*ms # Humphries, et al., 2006 
D2_WL_GPe_WL.w = np.random.choice([s,p,d],len(D2_WL_GPe_WL.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

############ GPe Projections 

# GPe SNR
GPeL_SNrL = Synapses(GPe_L,SNrL,weightEqs,on_pre=subW)
GPeL_SNrL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
GPeL_SNrL.delay = 3*ms # Humphries, et al., 2006 
GPeL_SNrL.w = np.random.choice([s,p],len(GPeL_SNrL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

GPeNL_SNrNL = Synapses(GPe_NL,SNrNL,weightEqs,on_pre=subW)
GPeNL_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
GPeNL_SNrNL.delay = 3*ms # Humphries, et al., 2006 
GPeNL_SNrNL.w = np.random.choice([s,p],len(GPeNL_SNrNL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

GPeWL_SNrWL = Synapses(GPe_WL,SNrWL,weightEqs,on_pre=subW)
GPeWL_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
GPeWL_SNrWL.delay = 3*ms # Humphries, et al., 2006 
GPeWL_SNrWL.w = np.random.choice([s,p],len(GPeWL_SNrWL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

# GPe to STN; Inhibitory
GPe_L_STN = Synapses(GPe_L,STN,weightEqs,subW)
GPe_L_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
GPe_L_STN.delay = 4*ms # Humphries, et al., 2006 
GPe_L_STN.w = np.random.choice([s,p,d],len(GPe_L_STN.i),p=[0.3,0.4,0.3]) # Humphries, et al., 2006  #rand(len(GPe_STN.i))

GPe_NL_STN = Synapses(GPe_NL,STN,weightEqs,subW)
GPe_NL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
GPe_NL_STN.delay = 4*ms # Humphries, et al., 2006 
GPe_NL_STN.w = np.random.choice([s,p,d],len(GPe_NL_STN.i),p=[0.3,0.4,0.3]) # Humphries, et al., 2006  #rand(len(GPe_STN.i))

GPe_WL_STN = Synapses(GPe_WL,STN,weightEqs,subW)
GPe_WL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
GPe_WL_STN.delay = 4*ms # Humphries, et al., 2006 
GPe_WL_STN.w = np.random.choice([s,p,d],len(GPe_WL_STN.i),p=[0.3,0.4,0.3]) # Humphries, et al., 2006  #rand(len(GPe_STN.i))

############ Poisson Projections 
# Poisson Cortex... introduce some noise into the whole system 
P_CortexL = Synapses(CortexL_Poisson,CortexL,weightEqs,on_pre=addW)
P_CortexL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
P_CortexL.delay = 5*ms
P_CortexL.w = 5

P_CortexNL = Synapses(CortexNL_Poisson,CortexNL,weightEqs,on_pre=addW)
P_CortexNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
P_CortexNL.delay = 5*ms
P_CortexNL.w = 5

P_CortexWL = Synapses(CortexWL_Poisson,CortexWL,weightEqs,on_pre=addW)
P_CortexWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
P_CortexWL.delay = 5*ms
P_CortexWL.w = 5

############ SNr Projections 
# Action
SNrL_Action =  Synapses(SNrL,LeverPress,on_pre='v-=10')
SNrL_Action.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
SNrL_Action.delay = 5*ms

SNrNL_NoAction = Synapses(SNrNL,NoLeverPress,on_pre='v-=10')
SNrNL_NoAction.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
SNrNL_NoAction.delay = 5*ms

SNrWL_WrongAction = Synapses(SNrWL,WrongLeverPress,on_pre='v-=10')
SNrWL_WrongAction.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
SNrWL_WrongAction.delay = 5*ms

# Thalamus
SNrL_ThalamusL = Synapses(SNrL,ThalamusL,weightEqs,on_pre=subW)
SNrL_ThalamusL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
SNrL_ThalamusL.delay = 5*ms#connect(j='k for k in range(i-w2, i+w2) if rand()<0.3', skip_if_invalid=True)
SNrL_ThalamusL.w = d

SNrNL_ThalamusNL = Synapses(SNrNL,ThalamusNL,weightEqs,on_pre=subW)
SNrNL_ThalamusNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
SNrNL_ThalamusNL.delay = 5*ms
SNrNL_ThalamusNL.w = d

SNrWL_ThalamusWL = Synapses(SNrWL,ThalamusWL,weightEqs,on_pre=subW)
SNrWL_ThalamusWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
SNrWL_ThalamusWL.delay = 5*ms
SNrWL_ThalamusWL.w = d

# Collaterals 
if  SNr_collaterals == 1:
    SNrL_SNrL = Synapses(SNrL,SNrL,weightEqs,on_pre=subW)
    SNrL_SNrL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrL_SNrL.delay = 1*ms # Humphries, et al., 2006
    SNrL_SNrL.w = np.random.choice([s,p],len(SNrL_SNrL.i),p=[0.5,0.5])
    
    SNrL_SNrNL = Synapses(SNrL,SNrNL,weightEqs,on_pre=subW)
    SNrL_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrL_SNrNL.delay = 1*ms # Humphries, et al., 2006
    SNrL_SNrNL.w = np.random.choice([s,p],len(SNrL_SNrNL.i),p=[0.5,0.5])
    
    SNrL_SNrWL = Synapses(SNrL,SNrWL,weightEqs,on_pre=subW)
    SNrL_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrL_SNrWL.delay = 1*ms # Humphries, et al., 2006
    SNrL_SNrWL.w = np.random.choice([s,p],len(SNrL_SNrWL.i),p=[0.5,0.5])
    
    SNrNL_SNrL = Synapses(SNrNL,SNrL,weightEqs,on_pre=subW)
    SNrNL_SNrL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrNL_SNrL.delay = 1*ms # Humphries, et al., 2006
    SNrNL_SNrL.w = np.random.choice([s,p],len(SNrNL_SNrL.i),p=[0.5,0.5])
    
    SNrNL_SNrNL = Synapses(SNrNL,SNrNL,weightEqs,on_pre=subW)
    SNrNL_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrNL_SNrNL.delay = 1*ms # Humphries, et al., 2006
    SNrNL_SNrNL.w = np.random.choice([s,p],len(SNrNL_SNrNL.i),p=[0.5,0.5])
    
    SNrNL_SNrWL = Synapses(SNrNL,SNrWL,weightEqs,on_pre=subW)
    SNrNL_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrNL_SNrWL.delay = 1*ms # Humphries, et al., 2006
    SNrNL_SNrWL.w = np.random.choice([s,p],len(SNrNL_SNrWL.i),p=[0.5,0.5])
    
    SNrWL_SNrL = Synapses(SNrWL,SNrL,weightEqs,on_pre=subW)
    SNrWL_SNrL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrWL_SNrL.delay = 1*ms # Humphries, et al., 2006
    SNrWL_SNrL.w = np.random.choice([s,p],len(SNrWL_SNrL.i),p=[0.5,0.5])
    
    SNrWL_SNrNL = Synapses(SNrWL,SNrNL,weightEqs,on_pre=subW)
    SNrWL_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrWL_SNrNL.delay = 1*ms # Humphries, et al., 2006
    SNrWL_SNrNL.w = np.random.choice([s,p],len(SNrWL_SNrNL.i),p=[0.5,0.5])
    
    SNrWL_SNrWL = Synapses(SNrWL,SNrWL,weightEqs,on_pre=subW)
    SNrWL_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
    SNrWL_SNrWL.delay = 1*ms # Humphries, et al., 2006
    SNrWL_SNrWL.w = np.random.choice([s,p],len(SNrWL_SNrWL.i),p=[0.5,0.5])

############ STN Projections 
# STN to SNr
STN_SNrL = Synapses(STN,SNrL,weightEqs,on_pre=addW)
STN_SNrL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
STN_SNrL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrL.w = p

STN_SNrNL = Synapses(STN,SNrNL,weightEqs,on_pre=addW)
STN_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
STN_SNrNL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrNL.w = p

STN_SNrWL = Synapses(STN,SNrWL,weightEqs,on_pre=addW)
STN_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
STN_SNrWL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrWL.w = p

# Excitatory STN to GPe
STN_GPe_L = Synapses(STN,GPe_L,weightEqs,on_pre=addW)
STN_GPe_L.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)  
STN_GPe_L.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_L.w = p # Humphries, et al., 2006... 

STN_GPe_NL = Synapses(STN,GPe_NL,weightEqs,on_pre=addW)
STN_GPe_NL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)  
STN_GPe_NL.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_NL.w = p # Humphries, et al., 2006... 

STN_GPe_WL = Synapses(STN,GPe_WL,weightEqs,on_pre=addW)
STN_GPe_WL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)  
STN_GPe_WL.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_WL.w = p # Humphries, et al., 2006... 

############ Thalamus Projections 

ThalCortexL = Synapses(ThalamusL,CortexL,weightEqs,on_pre=addW)
ThalCortexL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.3', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
ThalCortexL.delay = 5*ms 
ThalCortexL.w = 10*rand(len(ThalCortexL.i))

ThalCortexNL = Synapses(ThalamusNL,CortexNL,weightEqs,on_pre=addW)
ThalCortexNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.3', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
ThalCortexNL.delay = 5*ms 
ThalCortexNL.w = 10*rand(len(ThalCortexNL.i))

ThalCortexWL = Synapses(ThalamusWL,CortexWL,weightEqs,on_pre=addW)
ThalCortexWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.3', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-w2, i+w2) if rand()<0.33', skip_if_invalid=True)
ThalCortexWL.delay = 5*ms 
ThalCortexWL.w = 10*rand(len(ThalCortexWL.i))

#%% Record
ActionSpikes = SpikeMonitor(LeverPress)
NoActionSpikes = SpikeMonitor(NoLeverPress)
WrongActionSpikes = SpikeMonitor(WrongLeverPress)
InhSpikes = SpikeMonitor(InhInterL)
GPeSpikes = SpikeMonitor(GPe_L)
SNrLspikes = SpikeMonitor(SNrL)
SNrNLspikes = SpikeMonitor(SNrNL)
SNrWLspikes = SpikeMonitor(SNrWL)
STNspikes = SpikeMonitor(STN)
CortexLspikes= SpikeMonitor(CortexL)
CortexNLspikes = SpikeMonitor(CortexNL)
CortexWLspikes = SpikeMonitor(CortexWL)
DASpikes = SpikeMonitor(DA)
D2spikes = SpikeMonitor(D2_L)
D1_Lspikes = SpikeMonitor(D1_L)
D1_NLspikes = SpikeMonitor(D1_NL)
D1_WLspikes = SpikeMonitor(D1_WL)

if recordz == 1: 
    ActionTrace = StateMonitor(LeverPress,('v'),record=True)
    NoActionTrace = StateMonitor(NoLeverPress,('v'),record=True)
    WrongActionTrace = StateMonitor(NoLeverPress,('v'),record=True)
    
    D1_Lmon = StateMonitor(D1_L,('v'),record=True)
    GPeTrace = StateMonitor(GPe,('v'),record=True)
    SNrL_t = StateMonitor(SNrL,('v'),record=True)
    SNrNL_t = StateMonitor(SNrNL,('v'),record=True)

win = [0.25, 0.25, 0.5]

if learnAction == 1: 
   def calculate_FR(spike_monitor,integration_window,t): 
       bin_edges = np.arange((t-integration_window)/ms,t/ms,integration_window/4/ms)
       all_rates = []
       for i in range(0,len(bin_edges)-1):
           rate = np.size(np.where((spike_monitor.t > bin_edges[i]*ms)&(spike_monitor.t < bin_edges[i+1]*ms)))/n/integration_window
           all_rates.append(rate)
           smooth_rate = np.convolve(all_rates,win,mode='valid')
       return smooth_rate[len(smooth_rate)-1]*Hz
   

#   def DA_LTP(t,SpikeMon,SpikeMon2,SpikeMon3,SynapseMon,SynapseMon2):          
#       DAind = np.where(SpikeMon2.t > (t - window)) 
#       #SynapseMon.w = SynapseMon.w * 0.999
#       if len(DAind[0]) > 5: # was DA released?
#           CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
#           PresynapticInd = []
#           for i in range(0,len(np.unique(CortexIndex))):
#               pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
#               PresynapticInd.extend(pre[0])
#          act_MSN = np.where(SpikeMon3 > MSN_High) # find MSNs that are in a "high" state
#           high_synaptic_ind = np.concatenate([np.where(SynapseMon.j == x)[0] for x in act_MSN[0]]) # synaptic indicies projecting to high-state MSNs
#           s2 = set(high_synaptic_ind) # set object for high_synpatic ind
#          strengthen_synapse = [val for val in PresynapticInd if val in s2] # strengthen synapses that have MSNs in up state and cortical/DA input
#           not_strengthen_synapse = list(set(PresynapticInd) - set(strengthen_synapse))# weaken synapses that have glutamate but not upstate
#           SynapseMon.w[strengthen_synapse] +=   3*(SynapseMon.traceCon[strengthen_synapse] * mean(SynapseMon2.traceCon))     
#           SynapseMon.w[not_strengthen_synapse] -= (SynapseMon.w[not_strengthen_synapse] - s)/np.abs(SynapseMon.w[not_strengthen_synapse]/s)
#           SynapseMon.w = clip(SynapseMon.w, 0, wmax)
           
   def DA_LTP(t,SpikeMon,SpikeMon2,SpikeMon3,SynapseMon,SynapseMon2):          
       DAind = np.where(SpikeMon2.t > (t - window)) 
       #SynapseMon.w = SynapseMon.w * 0.999
       if len(DAind[0]) > 5: # was DA released?
           CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
           PresynapticInd = []
           for i in range(0,len(np.unique(CortexIndex))):
               pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
               PresynapticInd.extend(pre[0])
           act_MSN = np.where(SpikeMon3 > MSN_High) # find MSNs that are in a "high" state
           high_synaptic_ind = np.concatenate([np.where(SynapseMon.j == x)[0] for x in act_MSN[0]]) # synaptic indicies projecting to high-state MSNs
           s2 = set(high_synaptic_ind) # set object for high_synpatic ind
           strengthen_synapse = [val for val in PresynapticInd if val in s2] # strengthen synapses that have MSNs in up state and cortical/DA input
           not_strengthen_synapse = list(set(PresynapticInd) - set(strengthen_synapse))# weaken synapses that have glutamate but not upstate
           SynapseMon.w[strengthen_synapse] +=   1 
           SynapseMon.w[not_strengthen_synapse] -= (SynapseMon.w[not_strengthen_synapse] - s)/np.abs(SynapseMon.w[not_strengthen_synapse]/s)
           SynapseMon.w = clip(SynapseMon.w, 0, wmax)        

#network operations     
   rew_win=30*ms        
   @network_operation(dt=rew_win)
   def Reward(t):
       SNrL_rate = calculate_FR(SNrLspikes,integration_window,t)
       SNrNL_rate = calculate_FR(SNrNLspikes,integration_window,t)
       SNrWL_rate = calculate_FR(SNrWLspikes,integration_window,t)
       action = np.where(np.min([SNrL_rate,SNrNL_rate,SNrWL_rate]))
       if action[0] == 0:
          if SNrL_rate < SNr_thresh:
             DA.I += 2 # If a reward was recieved give DA
          else :
              DA.I = 0          
       else:
            DA.I = 0 # If no-reward was recieved zero out DA current 
           
   window = 20*ms
   @network_operation(dt=window)
   def calculate_LTP(t):
       DA_LTP(t,CortexLspikes,DASpikes,D1_Lspikes,CortexL_D1L,DA_D1L)    

   @network_operation(dt=window)
   def calculate_LTP2(t):
       DA_LTP(t,CortexNLspikes,DASpikes,D1_NLspikes,CortexNL_D1NL,DA_D1NL)           
    
   @network_operation(dt=window)    
   def calculate_LTP3(t):
       DA_LTP(t,CortexWLspikes,DASpikes,D1_WLspikes,CortexWL_D1WL,DA_D1WL)                   
   
#%% Run and analyze
       
if sequence == 1:  # reproduce figure 3 in humphries et al., 2006        
   # population activity
   ActionPop = PopulationRateMonitor(LeverPress)
   NoActionPop = PopulationRateMonitor(NoLeverPress)
   WrongActionPop = PopulationRateMonitor(WrongLeverPress)
   CortexPop = PopulationRateMonitor(CortexL)
   D1pop = PopulationRateMonitor(D1_L)
   GPePop = PopulationRateMonitor(GPe_L)
   SNrPop = PopulationRateMonitor(SNrL)
   SNrPopNL = PopulationRateMonitor(SNrNL)
   SNrPopWL = PopulationRateMonitor(SNrWL)
   STNpop = PopulationRateMonitor(STN)
   ThalamusPopL = PopulationRateMonitor(ThalamusL)
   ThalamusPopNL = PopulationRateMonitor(ThalamusNL)
   ThalamusPopWL = PopulationRateMonitor(ThalamusWL) 
   
   ThalamusL.I = 0
   ThalamusNL.I = 0
   ThalamusWL.I = 0
   DA.I = 0
   CortexL.I = 0
   run(sequence_duration,report='text')
   ThalamusL.I = 0
   ThalamusNL.I = 0
   ThalamusWL.I = 0
   DA.I = 0
   CortexL.I = 10
   CortexNL.I = 0
   run(sequence_duration,report='text')
   ThalamusL.I = 0
   ThalamusNL.I = 0
   ThalamusWL.I = 0
   DA.I = 0
   CortexL.I = 10
   CortexNL.I = 20
   run(sequence_duration,report='text')

   figure()
   plot(SNrPop.t/ms,SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(SNrPopNL.t/ms,SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(SNrPopWL.t/ms,SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('SNr Firing Rates')
   legend('R2U')
    
   figure()
   plot(ActionPop.t/ms,ActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(NoActionPop.t/ms,NoActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(WrongActionPop.t/ms,WrongActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('Action Firing Rates')
   legend('R2U')
   
   figure() 
   plot(D1pop.t/ms,D1pop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   xlabel('Time(ms)')
   ylabel('Firing Rate')


if popFiring == 1: # reproduce figure 2 in Humphries et al., 2006   
   # population activity
   GPePop = PopulationRateMonitor(GPe_L)
   SNrPop = PopulationRateMonitor(SNrL)
   SNrPopNL = PopulationRateMonitor(SNrNL)
   SNrPopWL = PopulationRateMonitor(SNrWL)
   STNpop = PopulationRateMonitor(STN)
   
   run (pop_duration,report='text') 

   STN_FR = mean(STNpop.smooth_rate(window='gaussian',width=binSize)/Hz)
   GPe_FR = mean(GPePop.smooth_rate(window='gaussian',width=binSize)/Hz)
   SNr_FR = mean(SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz)
   
   STNemp = 10;
   GPeEmp = 29;
   SNrEmp = 26

   print 'STN'   
   print STN_FR
   print 'GPe'
   print GPe_FR
   print 'SNr'
   print SNr_FR

   figure()
   #plt.legend('STN','GPe','SNr')
   plt.plot(1,STN_FR,'ob')
   plt.plot(2,GPe_FR,'ob')
   plt.plot(3,SNr_FR,'ob')
   plt.plot(1,STNemp,'g+',mew=10, ms=20)
   plt.plot(2,GPeEmp,'g+', mew=10, ms=20)
   plt.plot(3,SNrEmp,'g+', mew=10, ms=20)
   plt.errorbar(1,STNemp,(14-8)/2,capsize=11, elinewidth=2)
   plt.errorbar(2,GPeEmp,(36-26)/2,capsize=11, elinewidth=2)
   plt.errorbar(3,SNrEmp,(32-21)/2,capsize=11, elinewidth=2)
   xlim(0,4)
   ylim(8,40)
   title('Mean Tonic Firing Rates')
   
   print 'Population Firing Results'

if cortex_D1_action == 1:
    # population monitors 
    ActionPop = PopulationRateMonitor(LeverPress)
    NoActionPop = PopulationRateMonitor(NoLeverPress)
    WrongActionPop = PopulationRateMonitor(WrongLeverPress)
    CortexPop = PopulationRateMonitor(CortexL)
    D1pop = PopulationRateMonitor(D1_L)
    GPePop = PopulationRateMonitor(GPe_L)
    SNrPop = PopulationRateMonitor(SNrL)
    SNrPopNL = PopulationRateMonitor(SNrNL)
    SNrPopWL = PopulationRateMonitor(SNrWL)
    STNpop = PopulationRateMonitor(STN)
    ThalamusPopL = PopulationRateMonitor(ThalamusL)
    ThalamusPopNL = PopulationRateMonitor(ThalamusNL)
    ThalamusPopWL = PopulationRateMonitor(ThalamusWL)  
    
    ThalamusL.I = 0
    ThalamusNL.I = 0
    ThalamusWL.I = 0
    CortexL.I = 5
    CortexNL.I = 5
    CortexWL.I = 5
    run(cortex_D1_duration,report='text')
    
    # smoothed rates 
    SNrL_smoothFR = SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz
    SNrNL_smoothFR =SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz
    SNrWL_smoothFR =SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz
    
    # plots    
    figure()
    plot(ActionPop.t/ms,ActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    plot(NoActionPop.t/ms,NoActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
    plot(WrongActionPop.t/ms,WrongActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('Action Firing Rates')
    legend('RUU')
    
    figure()
    plot(SNrPop.t/ms,SNrL_smoothFR,'r')
    plot(SNrPopNL.t/ms,SNrNL_smoothFR,'b')
    plot(SNrPopWL.t/ms,SNrWL_smoothFR,'g')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('SNr Firing Rates')
    legend('RUU')
    
    figure()
    plot(CortexPop.t/ms,CortexPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('Cortical Firing Rates')
    legend('R')
    
    avg_D1L = np.mean(CortexL_D1L.w)
    avg_D1NL = np.mean(CortexNL_D1NL.w)
    avg_D1WL =  np.mean(CortexWL_D1WL.w)

    print 'Learned Action Results' 
    print np.str(np.sum(SNrPop.rate < action_thresh*Hz)) + '...rewarded action'
    print np.str(np.sum(SNrPopNL.rate < action_thresh*Hz)) + '...unrewared action'
    print np.str(np.sum(SNrPopWL.rate < action_thresh*Hz)) + '...unrewared action'

if learnAction == 1:
    # Population monitors     
    ActionPop = PopulationRateMonitor(LeverPress)
    NoActionPop = PopulationRateMonitor(NoLeverPress)
    WrongActionPop = PopulationRateMonitor(WrongLeverPress)
    CortexLpop = PopulationRateMonitor(CortexL)
    CortexNLpop = PopulationRateMonitor(CortexNL)
    CortexWLpop = PopulationRateMonitor(CortexWL)
    DApop = PopulationRateMonitor(DA)
    D1Lpop = PopulationRateMonitor(D1_L)
    D1NLpop = PopulationRateMonitor(D1_NL)
    D1WLpop = PopulationRateMonitor(D1_WL)
    GPePop = PopulationRateMonitor(GPe_L)
    SNrPop = PopulationRateMonitor(SNrL)
    SNrPopNL = PopulationRateMonitor(SNrNL)
    SNrPopWL = PopulationRateMonitor(SNrWL)
    STNpop = PopulationRateMonitor(STN)
    ThalamusPopL = PopulationRateMonitor(ThalamusL)
    ThalamusPopNL = PopulationRateMonitor(ThalamusNL)
    ThalamusPopWL = PopulationRateMonitor(ThalamusWL)  
    
    CortexL.I = 5
    CortexNL.I = 5
    CortexWL.I = 5
    run(learn_duration,report='text')
    
    figure()
    plot(ActionPop.t/ms,ActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    plot(NoActionPop.t/ms,NoActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
    plot(WrongActionPop.t/ms,WrongActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('Action Firing Rates')
    legend('R2U')
    
    figure()
    plot(SNrPop.t/ms,SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    plot(SNrPopNL.t/ms,SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
    plot(SNrPopWL.t/ms,SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('SNr Firing Rates')
    legend('R2U')
    
    print np.mean(CortexL_D1L.w)
    print np.mean(CortexNL_D1NL.w)
    print np.mean(CortexWL_D1WL.w)
    
    figure()
    plot(DApop.t/ms,DApop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('DA Firing Rates')
    
    print 'Learned Action Results' 
    print np.str(np.sum(SNrPop.rate < action_thresh*Hz)) + '...rewarded action'
    print np.str(np.sum(SNrPopNL.rate < action_thresh*Hz)) + '...unrewared action'
    print np.str(np.sum(SNrPopWL.rate < action_thresh*Hz)) + '...unrewared action'
    
    figure()
    plot(CortexLpop.t/ms,CortexLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    plot(CortexNLpop.t/ms,CortexNLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
    plot(CortexWLpop.t/ms,CortexWLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('Cortex Firing Rates')
    legend('R2U')
   
    figure()
    plot(D1Lpop.t/ms,D1Lpop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    plot(D1NLpop.t/ms,D1NLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
    plot(D1WLpop.t/ms,D1WLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('D1 Firing Rates')
    legend('R2U')   
   
if synfire == 1:
   Mgp = SpikeMonitor(CortexL) 
   CortexL.I = 0
   run(synfire_duration,report='text')
   
   figure()
   plot(Mgp.t/ms, 1.0*Mgp.i/group_size, '.')
   plot([0, synfire_duration/ms], np.arange(n_groups).repeat(2).reshape(-1, 2).T, 'k-')
   ylabel('group number')
   yticks(np.arange(n_groups))
   xlabel('time (ms)')
   title('Cortical Synfire for LevPress Channel')
   show()
   
   
if test_DA == 1:
   ActionPop = PopulationRateMonitor(LeverPress)
   NoActionPop = PopulationRateMonitor(NoLeverPress)
   WrongActionPop = PopulationRateMonitor(WrongLeverPress)
   CortexPop = PopulationRateMonitor(CortexL)
   DApop = PopulationRateMonitor(DA)
   D1pop = PopulationRateMonitor(D1_L)
   GPePop = PopulationRateMonitor(GPe_L)
   SNrPop = PopulationRateMonitor(SNrL)
   SNrPopNL = PopulationRateMonitor(SNrNL)
   SNrPopWL = PopulationRateMonitor(SNrWL)
   STNpop = PopulationRateMonitor(STN)
   ThalamusPopL = PopulationRateMonitor(ThalamusL)
   ThalamusPopNL = PopulationRateMonitor(ThalamusNL)
   ThalamusPopWL = PopulationRateMonitor(ThalamusWL)  

   DAvolts = StateMonitor(DA,('v'),record=True) 
    
   DA.I = 5
   run(DA_duration,report='text')
   DA.I = 0
   run(DA_duration,report='text')
   DA.I = 5
   run(DA_duration,report='text')

   figure()
   plot(DAvolts.t,DAvolts.v[0])    
    
   figure()
   plot(ActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(NoActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(WrongActionPop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('Action Firing Rates')
   legend('R2U')

   figure()
   plot(SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('SNr Firing Rates')
   legend('R2U')

   print  np.mean(DApop.rate)

   figure()
   plot(DApop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('DA Firing Rates')


