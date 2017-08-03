# -*- coding: utf-8 -*-
"""
Created on Wed Jul 05 11:49:34 2017

@author: elijah
"""


# removed STN scale
# got rid of +10 SNr neurons
# check list comprehension syntax with n.... seems odd that it goes from i-n to i+n (look at website)
# network should scale

from __future__ import division
from brian2 import *
import numpy as np
prefs.codegen.target = 'cython'

start_scope()     

# options
synfire = 0
SNr_collaterals = 1

# variables 
report_time = 60*second # how often to report simulation status 
pop_duration = 11000*ms # the duration to run simulations for population firing rates. This was 11 seconds in Humphries et al., 2006; 
sequence_duration = 1500*ms # As there are three stages this will result in a 3 seconds simulation
learn_duration = 10000*ms 
synfire_duration = 100*ms # a quick test to make sure the synfire chain is functioning correctly 
cortex_D1_duration = 3000*ms # a test of whether or not I can achieve more actions just through cortical-D1 plasticity 
DA_duration = 1000*ms # test to see if DA neurons fire
MSN_SNr_duration = 1000*ms # test to see how MSNs affect SNr activity
binSize = 50*ms 

# saving
save_root = 'C:/Users/elijah/Documents/model/BasalGanglia/FiguresFromTest/'

#%%
# Parameters
n = 100 # number of neurons

# TC Neuron
c0 = -63
b0 = 0.25
a0 = 0.02
d0 = 0.05

# RS Iz
c = -66
b= 0.2 
a = 0.02
d = 8

# FS Iz
c2 = -66
b2= 0.2
a2 = 0.1
d2 = 2  

# cortical
taum = 5*ms
taupsp = 0.325*ms#0.325*ms
weight = 4.86*mV

# STDP variables
taupre = taupost = 20*ms
wmax = 1500 # weights can go from 0 to 1
Apre = 0.01
Apost = -Apre*taupre/taupost*1.05

# Alpha synapses
alpha = 0
#tau_m=6*ms
#tau_ampa=20*ms
#g_synpk=0.01
#transwdth=1.0*ms
#memc = 2  # Membrane capacitance

# MSN plasticity
traceTau = traceTauPost = 20*ms  # this will change integration window # was 1200
traceConstant = 0.01 # this will change maximum weight change
tracePostConstant = -traceConstant*traceTau/traceTauPost*1.05

# Starting Weights... Different weights based upon synaptic distributions Humphries, et al., 2006 
s = 4 # Somatic 
p = 3 # proximal
d = 2 # distal 

# firing rate window
integration_window = 400*ms

# Cortical-MSN plasticity
MSN_High = -60 # mean value of up state (Plenz & kitai., 1998)
SNr_thresh = 35*Hz 

#lindhal et al., 2013 suggests the connectivity should be the same 
CortexD1_start = s

striatum_scale = 2 # how much to scale MSN numbers to GP/SNr inputs. This will account for convergance

#%% Neuron specific parameters

# GPe
GPe_spike = 30 # lindahl et al., 2013
GPe_a = 0.1 # Mandali et al., 2015
GPe_b = 0.2 # Mandali et al., 2015
GPe_c =-65 # lindahl et al., 2013
GPe_d = 2 # Mandali et al., 2015

# SNr 
SNr_spike = 30 # lindahl et al., 2013
SNr_a = 0.1 # Mandali et al., 2015
SNr_b = 0.2 # Mandali et al., 2015
SNr_c =-65 # Mandali et al., 2015 & # lindahl et al., 2013
SNr_d = 2 # Mandali et al., 2015

# STN
STN_spike = 30 # lindahl et al., 2013
STN_a = 0.005 # Mandali et al., 2015
STN_b = 0.265 # Mandali et al., 2015
STN_c =-65 # lindahl et al., 2013
STN_d = 1.5 # Mandali et al., 2015


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
du/dt = GPe_a*(GPe_b*v-u)/ms : 1 
I : 1
''' 
    
MSNeqs = ''' # Izikevich 2007 book 
dv/dt = (((v+80)*(v+25)-u+I)/ms)/50 : 1
du/dt = 0.01*(-20*(v+80)-u)/ms : 1
I : 1
'''

STNeqs = ''' # Mandali et al., 2015
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = STN_a*(STN_b*v-u)/ms : 1 
I : 1
'''  
 
if alpha == 1:
    SNrEqs = ''' # Mandali et al., 2015
    dv/dt = (((0.04)*v**2+(5)*v+140-u+I)/tau_m) - g_ampa/ms/memc : 1
    du/dt = SNr_a*(SNr_b*v-u)/ms : 1 
    dg_ampa/dt = -g_ampa/tau_ampa +z/ms : 1
    dz/dt = -z/tau_ampa +g_syn*Tr/ms : 1
    g_syn=g_synpk/(tau_ampa/ms*exp(-1)) : 1
    Tr=.25*(tanh((t/ms-tspike/ms)/.005)-tanh((t/ms-(tspike/ms +transwdth/ms))/.005)) : 1
    tspike : second
    I : 1
    '''   
else:
    SNrEqs = ''' # Mandali et al., 2015
    dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
    du/dt = SNr_a*(SNr_b*v-u)/ms : 1 
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
v=GPe_c
u+=GPe_d
'''

MSNreset = '''
v = -55 # 75 is mean value of down state (Plenz & kitai., 1998)
u+=150
'''

STNreset='''
v= STN_c
u+=STN_d
'''

if alpha == 1: 
   SNrReset='''
   v=SNr_c
   u+=SNr_d
   g_ampa=0
   tspike=t
   '''
else:
    SNrReset='''
    v=SNr_c
    u+=SNr_d
    '''

############ Plasticity Equations 
# weight equations
weightEqs='''
                w : 1
'''

subW='''
             v -= w
             '''

if alpha == 1:             
   D1_SNr='''
                 z_post -= g_synpk
                 '''             
else:             
    D1_SNr='''
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
             v_post += w
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

############ Cortical Neurons 
if synfire != 1:
   CortexL = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset0,method='euler')
    
   CortexL.v = c0
   CortexL.u = b0*c0
   CortexL.I = 0
    
   CortexNL = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset0,method='euler')
    
   CortexNL.v = c0
   CortexNL.u = b0*c0
   CortexNL.I = 0 

   CortexWL = NeuronGroup(n,model=eqs,threshold='v>30',reset=reset0,method='euler')

   CortexWL.v = c0
   CortexWL.u = b0*c0
   CortexWL.I = 0
   
elif synfire == 1:
     CortexL = NeuronGroup(N=n_groups*group_size, model=synfire_eqs,
                           threshold='v>30', reset=reset0, refractory=2*ms,
                           method='euler')
     CortexL.v = c0
     CortexL.u = b0*c0      
     CortexL.I = 0               
                             
     CortexNL = NeuronGroup(N=n_groups*group_size, model=synfire_eqs,
                           threshold='v>30', reset=reset0, refractory=2*ms,
                           method='euler')
     CortexNL.v = c0
     CortexNL.u = b0*c0     
     CortexNL.I = 0                           

     CortexWL = NeuronGroup(N=n_groups*group_size, model=synfire_eqs,
                            threshold='v>30', reset=reset0, refractory=2*ms,
                            method='euler')

     CortexWL.v = c0
     CortexWL.u = b0*c0       
     CortexWL.I = 0

     # input to the cortex 
     Pinput = SpikeGeneratorGroup(85, np.arange(85),
                                 np.random.randn(85)*1*ms + 50*ms)                       
                             
############ Striatal Neurons  
# D1 
D1_L = NeuronGroup(n*striatum_scale,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D1_L.v = c
D1_L.u = b*c
D1_L.I = 0

D1_NL = NeuronGroup(n*striatum_scale,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D1_NL.v = c
D1_NL.u = b*c
D1_NL.I = 0

D1_WL = NeuronGroup(n*striatum_scale,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D1_WL.v = c
D1_WL.u = b*c
D1_WL.I = 0

# D2
D2_L = NeuronGroup(n*striatum_scale,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D2_L.v = c
D2_L.u = b*c 
D2_L.I = 0 

D2_NL = NeuronGroup(n*striatum_scale,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D2_NL.v = c
D2_NL.u = b*c 
D2_NL.I = 0

D2_WL = NeuronGroup(n*striatum_scale,MSNeqs,threshold='v>40',reset=MSNreset,method='euler')

D2_WL.v = c
D2_WL.u = b*c 
D2_WL.I = 0

############ Dopamine Neurons 
DA = NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')
DA.v = c
DA.u = b*c
DA.I = 0

############ GPe Neurons
GPe_L = NeuronGroup(n,GPeEqs,threshold='v>GPe_spike',reset=GPeReset,method='euler') # Humphries, et al., 2006 

GPe_L.v = GPe_c
GPe_L.u = GPe_b*GPe_c
GPe_L.I = 0

GPe_NL = NeuronGroup(n,GPeEqs,threshold='v>GPe_spike',reset=GPeReset,method='euler') # Humphries, et al., 2006 

GPe_NL.v = GPe_c
GPe_NL.u = GPe_b*GPe_c
GPe_NL.I = 0

GPe_WL = NeuronGroup(n,GPeEqs,threshold='v>GPe_spike',reset=GPeReset,method='euler') # Humphries, et al., 2006 

GPe_WL.v = GPe_c
GPe_WL.u = GPe_b*GPe_c
GPe_WL.I = 0

############ SNr Neurons 
SNrL = NeuronGroup(n,SNrEqs,threshold='v>SNr_spike',reset=SNrReset,method='euler') # Humphries, et al., 2006 

SNrL.v = SNr_c
SNrL.u = SNr_b*SNr_c
SNrL.I = 0

SNrNL = NeuronGroup(n,SNrEqs,threshold='v>SNr_spike',reset=SNrReset,method='euler') # Humphries, et al., 2006 

SNrNL.v = SNr_c
SNrNL.u = SNr_b*SNr_c
SNrNL.I = 0

SNrWL = NeuronGroup(n,SNrEqs,threshold='v>SNr_spike',reset=SNrReset,method='euler') # Humphries, et al., 2006 

SNrWL.v = SNr_c
SNrWL.u = SNr_b*SNr_c
SNrWL.I = 0

############ STN Neurons 
STN = NeuronGroup(n,STNeqs,threshold='v>STN_spike',reset=STNreset,method='euler') # Humphries, et al., 2006 

STN.v = STN_c
STN.u = STN_c*STN_d
STN.I = 0

############ Thalamic Neurons 
ThalamusL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

ThalamusL.v = c0
ThalamusL.u = b0*c0
ThalamusL.I = 0

ThalamusNL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

ThalamusNL.v = c0
ThalamusNL.u = b0*c0
ThalamusNL.I = 0

ThalamusWL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

ThalamusWL.v = c0
ThalamusWL.u = b0*c0
ThalamusWL.I = 0

#%% Synapses

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

Cortex_MSN_prob = 30/n
CortexL_D1L = Synapses(CortexL,D1_L,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexL_D1L.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_MSN_prob', skip_if_invalid=True)
CortexL_D1L.delay = 10*ms # Humphries, et al., 2006 
CortexLD1start = CortexD1_start#*rand(len(CortexL_D1L.i))
CortexL_D1L.w = CortexD1_start #rand(len(Cortex_D1.i))

CortexNL_D1NL = Synapses(CortexNL,D1_NL,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexNL_D1NL.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_MSN_prob', skip_if_invalid=True)
CortexNL_D1NL.delay = 10*ms # Humphries, et al., 2006 
CortexNL_D1NL.w = CortexD1_start#s*rand(len(CortexNL_D1NL.i))

CortexWL_D1WL = Synapses(CortexWL,D1_WL,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexWL_D1WL.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_MSN_prob', skip_if_invalid=True)
CortexWL_D1WL.delay = 10*ms # Humphries, et al., 2006 
CortexWL_D1WL.w = CortexD1_start #s*rand(len(CortexWL_D1WL.i))

# Cortex D2 
CortexL_D2 = Synapses(CortexL,D2_L,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexL_D2.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_MSN_prob', skip_if_invalid=True)
CortexL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexL_D2.w = s

CortexNL_D2 = Synapses(CortexNL,D2_NL,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexNL_D2.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_MSN_prob', skip_if_invalid=True)
CortexNL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexNL_D2.w = s

CortexWL_D2 = Synapses(CortexWL,D2_WL,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexWL_D2.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_MSN_prob', skip_if_invalid=True)
CortexWL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexWL_D2.w = s

# Cortex STN - Hyperdirect pathway 
Cortex_STN_prob = 12.5/n
CortexL_STN = Synapses(CortexL,STN,weightEqs,on_pre=addW)
CortexL_STN.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_STN_prob', skip_if_invalid=True)
CortexL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexL_STN.w = d #d was suggested by humphries, but not enough activation 

CortexNL_STN = Synapses(CortexNL,STN,weightEqs,on_pre=addW)
CortexNL_STN.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_STN_prob', skip_if_invalid=True)
CortexNL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexNL_STN.w = d

CortexWL_STN = Synapses(CortexWL,STN,weightEqs,on_pre=addW)
CortexWL_STN.connect(j='k for k in range(i-n, i+n) if rand()<Cortex_STN_prob', skip_if_invalid=True)
CortexWL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexWL_STN.w = d

############ DopaminergicProjections 
# DA projections to D1
DA_MSN_prob = 30/n
DA_D1L = Synapses(DA,D1_L,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1L.connect(j='k for k in range(i-n, i+n) if rand()<DA_MSN_prob', skip_if_invalid=True) 
DA_D1L.delay = 5*ms
DA_D1L.w = d#rand(len(DA_D1L.i))

DA_D1NL = Synapses(DA,D1_NL,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1NL.connect(j='k for k in range(i-n, i+n) if rand()<DA_MSN_prob', skip_if_invalid=True) 
DA_D1NL.delay = 5*ms
DA_D1NL.w = d#rand(len(DA_D1NL.i))

DA_D1WL = Synapses(DA,D1_WL,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1WL.connect(j='k for k in range(i-n, i+n) if rand()<DA_MSN_prob', skip_if_invalid=True) 
DA_D1WL.delay = 5*ms
DA_D1WL.w = d#rand(len(DA_D1WL.i))

# DA projections to D1
DA_D2L = Synapses(DA,D2_L,DAstdp,on_pre=onpreDA_D2,on_post=onpostDA)
DA_D2L.connect(j='k for k in range(i-n, i+n) if rand()<DA_MSN_prob', skip_if_invalid=True) 
DA_D2L.delay = 5*ms
DA_D2L.w = d#rand(len(DA_D1L.i))

DA_D2NL = Synapses(DA,D2_NL,DAstdp,on_pre=onpreDA_D2,on_post=onpostDA)
DA_D2NL.connect(j='k for k in range(i-n, i+n) if rand()<DA_MSN_prob', skip_if_invalid=True) 
DA_D2NL.delay = 5*ms
DA_D2NL.w = d#rand(len(DA_D1NL.i))

DA_D2WL = Synapses(DA,D2_WL,DAstdp,on_pre=onpreDA_D2,on_post=onpostDA)
DA_D2WL.connect(j='k for k in range(i-n, i+n) if rand()<DA_MSN_prob', skip_if_invalid=True) 
DA_D2WL.delay = 5*ms
DA_D2WL.w = d#rand(len(DA_D1WL.i))

############ Striatal Projections 
# MSN-GP
MSN_GP = 30/n # this is the connectivity between D1-SNr and D2-GPe 
# D1 
D1_SNrL = Synapses(D1_L,SNrL,weightEqs,on_pre=D1_SNr)
D1_SNrL.connect(j='k for k in range(i-n, i+n) if rand()<MSN_GP', skip_if_invalid=True)
D1_SNrL.delay = 4*ms # Humphries, et al., 2006 
D1_SNrL.w = s# Humphries, et al., 2006  #rand(len(D1_SNrL.i))

D1_SNrNL = Synapses(D1_NL,SNrNL,weightEqs,on_pre=D1_SNr)
D1_SNrNL.connect(j='k for k in range(i-n, i+n) if rand()<MSN_GP', skip_if_invalid=True)
D1_SNrNL.delay = 4*ms # Humphries, et al., 2006 
D1_SNrNL.w = s #* convergance_scale # Humphries, et al., 2006  #rand(len(D1_SNrNL.i))

D1_SNrWL = Synapses(D1_WL,SNrWL,weightEqs,on_pre=D1_SNr)
D1_SNrWL.connect(j='k for k in range(i-n, i+n) if rand()<MSN_GP', skip_if_invalid=True)
D1_SNrWL.delay = 4*ms # Humphries, et al., 2006 
D1_SNrWL.w = s #* convergance_scale # Humphries, et al., 2006  #rand(len(D1_SNrWL.i))

# MSN collaterals
D1_collateral_prob = 10/100
D1L_NL = Synapses(D1_L,D1_NL,weightEqs,on_pre=subW)
D1L_NL.connect(j='k for k in range(i-n, i+n) if rand()<D1_collateral_prob', skip_if_invalid=True)
D1L_NL.delay = 2*ms # Humphries, et al., 2006 
D1L_NL.w = s# Humphries, et al., 2006  #rand(len(D1_SNrL.i))

D1L_WL = Synapses(D1_L,D1_WL,weightEqs,on_pre=subW)
D1L_WL.connect(j='k for k in range(i-n, i+n) if rand()<D1_collateral_prob', skip_if_invalid=True)
D1L_WL.delay = 2*ms # Humphries, et al., 2006 
D1L_WL.w = s# Humphries, et al., 2006  #rand(len(D1_SNrL.i))

D1NL_L = Synapses(D1_NL,D1_L,weightEqs,on_pre=subW)
D1NL_L.connect(j='k for k in range(i-n, i+n) if rand()<D1_collateral_prob', skip_if_invalid=True)
D1NL_L.delay = 2*ms # Humphries, et al., 2006 
D1NL_L.w = s# Humphries, et al., 2006  #rand(len(D1_SNrL.i))

D1NL_WL = Synapses(D1_NL,D1_WL,weightEqs,on_pre=subW)
D1NL_WL.connect(j='k for k in range(i-n, i+n) if rand()<D1_collateral_prob', skip_if_invalid=True)
D1NL_WL.delay = 2*ms # Humphries, et al., 2006 
D1NL_WL.w = s# Humphries, et al., 2006  #rand(len(D1_SNrL.i))

D1WL_L = Synapses(D1_WL,D1_L,weightEqs,on_pre=subW)
D1WL_L.connect(j='k for k in range(i-n, i+n) if rand()<D1_collateral_prob', skip_if_invalid=True)
D1WL_L.delay = 2*ms # Humphries, et al., 2006 
D1WL_L.w = s# Humphries, et al., 2006  #rand(len(D1_SNrL.i))

D1WL_NL = Synapses(D1_WL,D1_NL,weightEqs,on_pre=subW)
D1WL_NL.connect(j='k for k in range(i-n, i+n) if rand()<D1_collateral_prob', skip_if_invalid=True)
D1WL_NL.delay = 2*ms # Humphries, et al., 2006 
D1WL_NL.w = s# Humphries, et al., 2006  #rand(len(D1_SNrL.i))

# D2 
D2_L_GPe_L = Synapses(D2_L,GPe_L,weightEqs,on_pre=subW)
D2_L_GPe_L.connect(j='k for k in range(i-n, i+n) if rand()<MSN_GP', skip_if_invalid=True) # Lindhal et al., 2013 suggest keeping D1-SNr and D2-GPe the same 
D2_L_GPe_L.delay = 5*ms # Humphries, et al., 2006 
D2_L_GPe_L.w = np.random.choice([s,p,d],len(D2_L_GPe_L.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

D2_NL_GPe_NL = Synapses(D2_NL,GPe_NL,weightEqs,on_pre=subW)
D2_NL_GPe_NL.connect(j='k for k in range(i-n, i+n) if rand()<MSN_GP', skip_if_invalid=True)
D2_NL_GPe_NL.delay = 5*ms # Humphries, et al., 2006 
D2_NL_GPe_NL.w = np.random.choice([s,p,d],len(D2_NL_GPe_NL.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

D2_WL_GPe_WL = Synapses(D2_WL,GPe_WL,weightEqs,on_pre=subW)
D2_WL_GPe_WL.connect(j='k for k in range(i-n, i+n) if rand()<MSN_GP', skip_if_invalid=True)
D2_WL_GPe_WL.delay = 5*ms # Humphries, et al., 2006 
D2_WL_GPe_WL.w = np.random.choice([s,p,d],len(D2_WL_GPe_WL.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

############ GPe Projections 

# GPe SNR
GPe_SNr_prob = 25/n
GPeL_SNrL = Synapses(GPe_L,SNrL,weightEqs,on_pre=subW)
GPeL_SNrL.connect(j='k for k in range(i-n, i+n) if rand()<GPe_SNr_prob', skip_if_invalid=True) 
GPeL_SNrL.delay = 3*ms # Humphries, et al., 2006 
GPeL_SNrL.w = np.random.choice([s,p],len(GPeL_SNrL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

GPeNL_SNrNL = Synapses(GPe_NL,SNrNL,weightEqs,on_pre=subW)
GPeNL_SNrNL.connect(j='k for k in range(i-n, i+n) if rand()<GPe_SNr_prob', skip_if_invalid=True) 
GPeNL_SNrNL.delay = 3*ms # Humphries, et al., 2006 
GPeNL_SNrNL.w = np.random.choice([s,p],len(GPeNL_SNrNL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

GPeWL_SNrWL = Synapses(GPe_WL,SNrWL,weightEqs,on_pre=subW)
GPeWL_SNrWL.connect(j='k for k in range(i-n, i+n) if rand()<GPe_SNr_prob', skip_if_invalid=True) 
GPeWL_SNrWL.delay = 3*ms # Humphries, et al., 2006 
GPeWL_SNrWL.w = np.random.choice([s,p],len(GPeWL_SNrWL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

# GPe to STN; Inhibitory
GPe_STN_prob = 25/n
GPe_L_STN = Synapses(GPe_L,STN,weightEqs,subW)
GPe_L_STN.connect(j='k for k in range(i-n, i+n) if rand()<GPe_STN_prob', skip_if_invalid=True) 
GPe_L_STN.delay = 5*ms # Lindahl et al., 2013 
GPe_L_STN.w = np.random.choice([s,p,d],len(GPe_L_STN.i),p=[0.3,0.4,0.3]) # Humphries, et al., 2006  #rand(len(GPe_STN.i))

GPe_NL_STN = Synapses(GPe_NL,STN,weightEqs,subW)
GPe_NL_STN.connect(j='k for k in range(i-n, i+n) if rand()<GPe_STN_prob', skip_if_invalid=True) 
GPe_NL_STN.delay = 5*ms # Lindahl et al., 2013 
GPe_NL_STN.w = np.random.choice([s,p,d],len(GPe_NL_STN.i),p=[0.3,0.4,0.3]) # Humphries, et al., 2006  #rand(len(GPe_STN.i))

GPe_WL_STN = Synapses(GPe_WL,STN,weightEqs,subW)
GPe_WL_STN.connect(j='k for k in range(i-n, i+n) if rand()<GPe_STN_prob', skip_if_invalid=True) 
GPe_WL_STN.delay = 5*ms # Lindahl et al., 2013 
GPe_WL_STN.w = np.random.choice([s,p,d],len(GPe_WL_STN.i),p=[0.3,0.4,0.3]) # Humphries, et al., 2006  #rand(len(GPe_STN.i))

############ SNr Projections 

# Thalamus
SNr_Thalamus_prob = 12.5/n
SNrL_ThalamusL = Synapses(SNrL,ThalamusL,weightEqs,on_pre=subW)
SNrL_ThalamusL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_Thalamus_prob', skip_if_invalid=True)
SNrL_ThalamusL.delay = 5*ms
SNrL_ThalamusL.w = d

SNrNL_ThalamusNL = Synapses(SNrNL,ThalamusNL,weightEqs,on_pre=subW)
SNrNL_ThalamusNL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_Thalamus_prob', skip_if_invalid=True)
SNrNL_ThalamusNL.delay = 5*ms
SNrNL_ThalamusNL.w = d

SNrWL_ThalamusWL = Synapses(SNrWL,ThalamusWL,weightEqs,on_pre=subW)
SNrWL_ThalamusWL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_Thalamus_prob', skip_if_invalid=True)
SNrWL_ThalamusWL.delay = 5*ms
SNrWL_ThalamusWL.w = d

# Collaterals 
if  SNr_collaterals == 1:
    SNr_prob = 25/n # Humphries et al., 2006 uses 0.25/n... n is less than 100 though 
    SNrL_SNrL = Synapses(SNrL,SNrL,weightEqs,on_pre=subW)
    SNrL_SNrL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrL_SNrL.delay = 1*ms # Humphries, et al., 2006
    SNrL_SNrL.w = np.random.choice([s,p],len(SNrL_SNrL.i),p=[0.5,0.5])
    
    SNrL_SNrNL = Synapses(SNrL,SNrNL,weightEqs,on_pre=subW)
    SNrL_SNrNL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrL_SNrNL.delay = 1*ms # Humphries, et al., 2006
    SNrL_SNrNL.w = np.random.choice([s,p],len(SNrL_SNrNL.i),p=[0.5,0.5])
    
    SNrL_SNrWL = Synapses(SNrL,SNrWL,weightEqs,on_pre=subW)
    SNrL_SNrWL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrL_SNrWL.delay = 1*ms # Humphries, et al., 2006
    SNrL_SNrWL.w = np.random.choice([s,p],len(SNrL_SNrWL.i),p=[0.5,0.5])
    
    SNrNL_SNrL = Synapses(SNrNL,SNrL,weightEqs,on_pre=subW)
    SNrNL_SNrL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrNL_SNrL.delay = 1*ms # Humphries, et al., 2006
    SNrNL_SNrL.w = np.random.choice([s,p],len(SNrNL_SNrL.i),p=[0.5,0.5])
    
    SNrNL_SNrNL = Synapses(SNrNL,SNrNL,weightEqs,on_pre=subW)
    SNrNL_SNrNL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrNL_SNrNL.delay = 1*ms # Humphries, et al., 2006
    SNrNL_SNrNL.w = np.random.choice([s,p],len(SNrNL_SNrNL.i),p=[0.5,0.5])
    
    SNrNL_SNrWL = Synapses(SNrNL,SNrWL,weightEqs,on_pre=subW)
    SNrNL_SNrWL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrNL_SNrWL.delay = 1*ms # Humphries, et al., 2006
    SNrNL_SNrWL.w = np.random.choice([s,p],len(SNrNL_SNrWL.i),p=[0.5,0.5])
    
    SNrWL_SNrL = Synapses(SNrWL,SNrL,weightEqs,on_pre=subW)
    SNrWL_SNrL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrWL_SNrL.delay = 1*ms # Humphries, et al., 2006
    SNrWL_SNrL.w = np.random.choice([s,p],len(SNrWL_SNrL.i),p=[0.5,0.5])
    
    SNrWL_SNrNL = Synapses(SNrWL,SNrNL,weightEqs,on_pre=subW)
    SNrWL_SNrNL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrWL_SNrNL.delay = 1*ms # Humphries, et al., 2006
    SNrWL_SNrNL.w = np.random.choice([s,p],len(SNrWL_SNrNL.i),p=[0.5,0.5])
    
    SNrWL_SNrWL = Synapses(SNrWL,SNrWL,weightEqs,on_pre=subW)
    SNrWL_SNrWL.connect(j='k for k in range(i-n, i+n) if rand()<SNr_prob', skip_if_invalid=True)
    SNrWL_SNrWL.delay = 1*ms # Humphries, et al., 2006
    SNrWL_SNrWL.w = np.random.choice([s,p],len(SNrWL_SNrWL.i),p=[0.5,0.5])

############ STN Projections 
# STN to SNr
STN_SNr_prob = 85/n
STN_SNrL = Synapses(STN,SNrL,weightEqs,on_pre=addW)
STN_SNrL.connect(j='k for k in range(i-n, i+n) if rand()<STN_SNr_prob', skip_if_invalid=True)
STN_SNrL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrL.w = p

STN_SNrNL = Synapses(STN,SNrNL,weightEqs,on_pre=addW)
STN_SNrNL.connect(j='k for k in range(i-n, i+n) if rand()<STN_SNr_prob', skip_if_invalid=True)
STN_SNrNL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrNL.w = p

STN_SNrWL = Synapses(STN,SNrWL,weightEqs,on_pre=addW)
STN_SNrWL.connect(j='k for k in range(i-n, i+n) if rand()<STN_SNr_prob', skip_if_invalid=True)
STN_SNrWL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrWL.w = p

# Excitatory STN to GPe
STN_GPe_prob = 60/n
STN_GPe_L = Synapses(STN,GPe_L,weightEqs,on_pre=addW)
STN_GPe_L.connect(j='k for k in range(i-n, i+n) if rand()<STN_GPe_prob', skip_if_invalid=True)  
STN_GPe_L.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_L.w = p # Humphries, et al., 2006... 

STN_GPe_NL = Synapses(STN,GPe_NL,weightEqs,on_pre=addW)
STN_GPe_NL.connect(j='k for k in range(i-n, i+n) if rand()<STN_GPe_prob', skip_if_invalid=True)  
STN_GPe_NL.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_NL.w = p # Humphries, et al., 2006... 

STN_GPe_WL = Synapses(STN,GPe_WL,weightEqs,on_pre=addW)
STN_GPe_WL.connect(j='k for k in range(i-n, i+n) if rand()<STN_GPe_prob', skip_if_invalid=True)  
STN_GPe_WL.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_WL.w = p # Humphries, et al., 2006... 

############ Thalamus Projections 
Thal_Cortex_prob = 12.5/n
ThalCortexL = Synapses(ThalamusL,CortexL,weightEqs,on_pre=addW)
ThalCortexL.connect(j='k for k in range(i-n, i+n) if rand()<Thal_Cortex_prob', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-n, i+n) if rand()<0.33', skip_if_invalid=True)
ThalCortexL.delay = 5*ms 
ThalCortexL.w = p

ThalCortexNL = Synapses(ThalamusNL,CortexNL,weightEqs,on_pre=addW)
ThalCortexNL.connect(j='k for k in range(i-n, i+n) if rand()<Thal_Cortex_prob', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-n, i+n) if rand()<0.33', skip_if_invalid=True)
ThalCortexNL.delay = 5*ms 
ThalCortexNL.w = p

ThalCortexWL = Synapses(ThalamusWL,CortexWL,weightEqs,on_pre=addW)
ThalCortexWL.connect(j='k for k in range(i-n, i+n) if rand()<Thal_Cortex_prob', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-n, i+n) if rand()<0.33', skip_if_invalid=True)
ThalCortexWL.delay = 5*ms 
ThalCortexWL.w = p

Thalamus_D1_prob = 12.5/n
Thal_D1L = Synapses(ThalamusL,D1_L,weightEqs,on_pre=addW)
Thal_D1L.connect(j='k for k in range(i-n, i+n) if rand()<Thalamus_D1_prob', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-n, i+n) if rand()<0.33', skip_if_invalid=True)
Thal_D1L.delay = 5*ms 
Thal_D1L.w = p

Thal_D1NL = Synapses(ThalamusL,D1_NL,weightEqs,on_pre=addW)
Thal_D1NL.connect(j='k for k in range(i-n, i+n) if rand()<Thalamus_D1_prob', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-n, i+n) if rand()<0.33', skip_if_invalid=True)
Thal_D1NL.delay = 5*ms 
Thal_D1NL.w = p

Thal_D1WL = Synapses(ThalamusL,D1_WL,weightEqs,on_pre=addW)
Thal_D1WL.connect(j='k for k in range(i-n, i+n) if rand()<Thalamus_D1_prob', skip_if_invalid=True) #connect(i=[0,1,2],j=[0,1,2])# j='k for k in range(i-n, i+n) if rand()<0.33', skip_if_invalid=True)
Thal_D1WL.delay = 5*ms 
Thal_D1WL.w = p