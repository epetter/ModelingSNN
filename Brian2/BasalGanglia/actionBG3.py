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
# EP 161214
##############

#%% import and setup 
from brian2 import *
import numpy as np
import scipy as sp
import scipy.stats
#prefs.codegen.target = 'weave'
start_scope()
duration =10000*ms #50000*ms#1000000*ms

recordz = 0
plotz = 0
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
traceTau = traceTauPost = 1200*ms  # this will change integration window
traceConstant = 0.01 # this will change maximum weight change
tracePostConstant = -traceConstant*traceTau/traceTauPost*1.05

# Starting Weights... Different weights based upon synaptic distributions Humphries, et al., 2006 
s = 10 # Somatic 
p = 8 # proximal
d = 5 # distal 

#%% Equations 
############ Neuron Eqs

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

############ WTA neurons 
CortexPoisson = PoissonGroup(n, np.arange(n)*Hz + 10*Hz) # input to STN neurons
#GPePoisson = PoissonGroup(n, np.arange(n)*Hz + 10*Hz) # input to STN neurons

############ WTA neurons 
LeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
NoLeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
WrongLeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
#InhInter = NeuronGroup(5,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterNL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterWL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')

############ Cortical Neurons 
CortexL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

CortexL.v = c0
CortexL.u = b0*c0

CortexNL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

CortexNL.v = c0
CortexNL.u = b0*c0

CortexWL = NeuronGroup(n,eqs,threshold='v>30',reset=reset0,method='euler')

CortexWL.v = c0
CortexWL.u = b0*c0

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

#STN_NL = NeuronGroup(n,eqs2,threshold='v>20',reset=reset,method='euler') # Humphries, et al., 2006 

#STN_NL.v = c
#STN_NL.u = b*c

#STN_WL = NeuronGroup(n,eqs2,threshold='v>20',reset=reset,method='euler') # Humphries, et al., 2006 

#STN_WL.v = c
#STN_WL.u = b*c

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
# Cortex WTA
CortexL_Lever = Synapses(CortexL,LeverPress,on_pre='v+=25')
CortexL_Lever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) #i=[0,1,2],j=0)#j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
CortexL_Lever.delay = 15*ms

CortexNL_NoLever = Synapses(CortexNL,NoLeverPress,on_pre='v+=25')
CortexNL_NoLever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) #(i=[3,4,5],j=0)#j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
CortexNL_NoLever.delay = 15*ms

CortexWL_WrongLever = Synapses(CortexWL,WrongLeverPress,on_pre='v+=25')
CortexWL_WrongLever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) #i=[6,7,8],j=0)#j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
CortexWL_WrongLever.delay = 15*ms

# Cortex D1
CortexL_D1L = Synapses(CortexL,D1_L,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexL_D1L.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
CortexL_D1L.delay = 10*ms # Humphries, et al., 2006 
CortexLD1start = 100#s#*rand(len(CortexL_D1L.i))
CortexL_D1L.w = CortexLD1start #rand(len(Cortex_D1.i))

CortexNL_D1NL = Synapses(CortexNL,D1_NL,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexNL_D1NL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
CortexNL_D1NL.delay = 10*ms # Humphries, et al., 2006 
CortexNL_D1NL.w = s#s*rand(len(CortexNL_D1NL.i))

CortexWL_D1WL = Synapses(CortexWL,D1_WL,MSNstdp,on_pre=MSNpre,on_post=MSNpost) #on_pre='v=+10')
CortexWL_D1WL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
CortexWL_D1WL.delay = 10*ms # Humphries, et al., 2006 
CortexWL_D1WL.w = s#s*rand(len(CortexWL_D1WL.i))

# Cortex D2 
CortexL_D2 = Synapses(CortexL,D2_L,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexL_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.6', skip_if_invalid=True)
CortexL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexL_D2.w = s

CortexNL_D2 = Synapses(CortexNL,D2_NL,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexNL_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.6', skip_if_invalid=True)
CortexNL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexNL_D2.w = s

CortexWL_D2 = Synapses(CortexWL,D2_WL,weightEqs,on_pre=addW) #on_pre='v=+10')
CortexWL_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.6', skip_if_invalid=True)
CortexWL_D2.delay = 10*ms # Humphries, et al., 2006 
CortexWL_D2.w = s

# Cortex STN - Hyperdirect pathway 
CortexL_STN = Synapses(CortexL,STN,weightEqs,on_pre=addW)
CortexL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)
CortexL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexL_STN.w = d #d was suggested by humphries, but not enough activation 

CortexNL_STN = Synapses(CortexNL,STN,weightEqs,on_pre=addW)
CortexNL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)
CortexNL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexNL_STN.w = d

CortexWL_STN = Synapses(CortexWL,STN,weightEqs,on_pre=addW)
CortexWL_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)
CortexWL_STN.delay = 2.5*ms # Humphries, et al., 2006 
CortexWL_STN.w = d

############ DopaminergicProjections 
# DA projections to D1
DA_D1L = Synapses(DA,D1_L,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1L.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
#DA_D1.connect(condition='i!=j',p=0.75,skip_if_invalid=True)  #p=0.1) # I think I'd rather have a broad range of DA with synapse specific plasticity at places with glut and DA
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
#DA_D1.connect(condition='i!=j',p=0.75,skip_if_invalid=True)  #p=0.1) # I think I'd rather have a broad range of DA with synapse specific plasticity at places with glut and DA
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
D2_L_GPe_L.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)
D2_L_GPe_L.delay = 5*ms # Humphries, et al., 2006 
D2_L_GPe_L.w = np.random.choice([s,p,d],len(D2_L_GPe_L.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

D2_NL_GPe_NL = Synapses(D2_NL,GPe_NL,weightEqs,on_pre=subW)
D2_NL_GPe_NL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)
D2_NL_GPe_NL.delay = 5*ms # Humphries, et al., 2006 
D2_NL_GPe_NL.w = np.random.choice([s,p,d],len(D2_NL_GPe_NL.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

D2_WL_GPe_WL = Synapses(D2_WL,GPe_WL,weightEqs,on_pre=subW)
D2_WL_GPe_WL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)
D2_WL_GPe_WL.delay = 5*ms # Humphries, et al., 2006 
D2_WL_GPe_WL.w = np.random.choice([s,p,d],len(D2_WL_GPe_WL.i),p=[0.33,0.33,0.34]) # Humphries, et al., 2006 

############ GPe Projections 
# Collaterals
#GPe_GPe = Synapses(GPe,GPe,weightEqs,on_pre=subW)
#GPe_GPe.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
#GPe_GPe.delay = 4*ms # Humphries, et al., 2006 
#GPe_GPe.w = np.random.choice([s,p]) # Humphries, et al., 2006 

# GPe SNR
GPeL_SNrL = Synapses(GPe_L,SNrL,weightEqs,on_pre=subW)
GPeL_SNrL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.45', skip_if_invalid=True) 
GPeL_SNrL.delay = 3*ms # Humphries, et al., 2006 
GPeL_SNrL.w = np.random.choice([s,p],len(GPeL_SNrL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

GPeNL_SNrNL = Synapses(GPe_NL,SNrNL,weightEqs,on_pre=subW)
GPeNL_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.45', skip_if_invalid=True) 
GPeNL_SNrNL.delay = 3*ms # Humphries, et al., 2006 
GPeNL_SNrNL.w = np.random.choice([s,p],len(GPeNL_SNrNL.i),p=[0.5,0.5]) # Humphries, et al., 2006 

GPeWL_SNrWL = Synapses(GPe_WL,SNrWL,weightEqs,on_pre=subW)
GPeWL_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.45', skip_if_invalid=True) 
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
#P_CortexL = Synapses(CortexPoisson,CortexL,weightEqs,on_pre=addW)
#P_CortexL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
#P_CortexL.delay = 5*ms
#P_CortexL.w = 5

#P_CortexNL = Synapses(CortexPoisson,CortexNL,weightEqs,on_pre=addW)
#P_CortexNL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
#P_CortexNL.delay = 5*ms
#P_CortexNL.w = 5

#P_CortexWL = Synapses(CortexPoisson,CortexWL,weightEqs,on_pre=addW)
#P_CortexWL.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
#P_CortexWL.delay = 5*ms
#P_CortexWL.w = 5

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
STN_SNrL.w = d

STN_SNrNL = Synapses(STN,SNrNL,weightEqs,on_pre=addW)
STN_SNrNL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
STN_SNrNL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrNL.w = d

STN_SNrWL = Synapses(STN,SNrWL,weightEqs,on_pre=addW)
STN_SNrWL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)
STN_SNrWL.delay = 1.5*ms # Humphries, et al., 2006 
STN_SNrWL.w = d

# Excitatory STN to GPe
STN_GPe_L = Synapses(STN,GPe_L,weightEqs,on_pre=addW)
STN_GPe_L.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)  
STN_GPe_L.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_L.w = d # Humphries, et al., 2006... added because GPe wasn't spiking 

STN_GPe_NL = Synapses(STN,GPe_NL,weightEqs,on_pre=addW)
STN_GPe_NL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)  
STN_GPe_NL.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_NL.w = d # Humphries, et al., 2006... added because GPe wasn't spiking 

STN_GPe_WL = Synapses(STN,GPe_WL,weightEqs,on_pre=addW)
STN_GPe_WL.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True)  
STN_GPe_WL.delay = 2*ms # Humphries, et al., 2006 
STN_GPe_WL.w = d # Humphries, et al., 2006... added because GPe wasn't spiking 

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
CortexLSpikes= SpikeMonitor(CortexL)
CortexNLspikes = SpikeMonitor(CortexNL)
CortexWLspikes = SpikeMonitor(CortexWL)
DASpikes = SpikeMonitor(DA)
D2spikes = SpikeMonitor(D2_L)
D1spikes = SpikeMonitor(D1_L)
D1nlSpikes = SpikeMonitor(D1_NL)
#GPeTrace = StateMonitor(GPe,('v'),record=True)

#CortexD1Lmon= StateMonitor(Cortex_D1L,('w'),record=True)
GPePop = PopulationRateMonitor(GPe_L)
SNrPop = PopulationRateMonitor(SNrL)
SNrPopNL = PopulationRateMonitor(SNrNL)
SNrPopWL = PopulationRateMonitor(SNrWL)
STNpop = PopulationRateMonitor(STN)

# population activity
ThalamusPopL = PopulationRateMonitor(ThalamusL)
ThalamusPopNL = PopulationRateMonitor(ThalamusNL)
ThalamusPopWL = PopulationRateMonitor(ThalamusWL)
CortexPop = PopulationRateMonitor(CortexL)
D1pop = PopulationRateMonitor(D1_L)

if recordz == 1: 
    ActionTrace = StateMonitor(LeverPress,('v'),record=True)
    NoActionTrace = StateMonitor(NoLeverPress,('v'),record=True)
    WrongActionTrace = StateMonitor(NoLeverPress,('v'),record=True)
    #Action2Trace = StateMonitor(LeverPress2,('v'),record=True)
    
    SNrL_t = StateMonitor(SNrL,('v'),record=True)
    SNrNL_t = StateMonitor(SNrNL,('v'),record=True)
    
    D1_Lmon = StateMonitor(D1_L,('v'),record=True)

#%% Network Operator
#@network_operation(dt=50*ms)
#def updateWeights():
#    actionSelection = np.array([len(ActionSpikes.t[ActionSpikes.t > defaultclock.t-50*ms]),len(NoActionSpikes.t[NoActionSpikes.t > defaultclock.t-50*ms]),len(WrongActionSpikes.t[WrongActionSpikes.t > defaultclock.t-50*ms])])      
#    action = np.where(actionSelection == np.max(actionSelection))
#    #print actionSelection
#    if len(action[0])<2:
#        if action[0]<1:
#            CortexL_D1L.w += 1
           #indz = np.where(Cortex_D1.t > (defaultclock.t - 50*ms)) #Cortex_D2_DAtrace.t > defaultclock.t-10*ms)
           #Cortex_D1.w[indz[0]] += 1
          #print 'learn'

@network_operation(dt=500*ms) # update the sensory input every 100ms
def sensInput():
    if ThalamusL.I[0] > 0: # if one is above zero all should be above zero 
       ThalamusL.I[:] = 0
       ThalamusNL.I[:] = 0
       ThalamusWL.I[:] = 0
    else: 
        ThalamusL.I[:] = 5
        ThalamusNL.I[:] = 5
        ThalamusWL.I[:] = 5

rewWin=150*ms        
@network_operation(dt=rewWin)
def Reward():
    actionSelection = np.array([len(ActionSpikes.t[ActionSpikes.t > defaultclock.t-rewWin]),len(NoActionSpikes.t[NoActionSpikes.t > defaultclock.t-50*ms]),len(WrongActionSpikes.t[WrongActionSpikes.t > defaultclock.t-50*ms])])      
    action = np.where(actionSelection == np.max(actionSelection))
    if len(action[0])<2:
        if action[0]<1:
            if ThalamusL.I[0] > 0:
                DA.v += 20
            #print 'here'

window = 50*ms
@network_operation(dt=window)
def DA_LTP():
    SpikeMon=CortexLSpikes
    SpikeMon2=DASpikes
    SynapseMon=CortexL_D1L
    SynapseMon2=DA_D1L
    DAind = np.where(SpikeMon2.t > (defaultclock.t - window)) 
    SynapseMon.w = SynapseMon.w * 0.999
    if len(DAind[0]) > 5: # was DA released?
        CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
        PresynapticInd = []
        for i in range(0,len(np.unique(CortexIndex))):
            pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
            PresynapticInd.extend(pre[0])
        PostsynapticInd = np.where(SynapseMon.j[PresynapticInd])  
        NotPostsynapticInd = np.delete(range(0,len(SynapseMon.j)),PostsynapticInd[0],0)
        #PostsynapticNeuron = np.unique(PostsynapticInd) # DA firing is broad, but need this if I want to figure out which neurons recieve both
        #DA ... # DA firing is broad, but need this if I want to figure out which neurons recieve both
        SynapseMon.w[PostsynapticInd[0]] += 3*(SynapseMon.traceCon[PostsynapticInd[0]] * mean(SynapseMon2.traceCon))        
        SynapseMon.w[NotPostsynapticInd[0]] += 3*(SynapseMon.traceConPost[NotPostsynapticInd[0]] * mean(SynapseMon2.traceCon))
        SynapseMon.w = clip(SynapseMon.w, 0, wmax)
 
@network_operation(dt=window)
def DA_LTP2():
    SpikeMon=CortexNLspikes
    SpikeMon2=DASpikes
    SynapseMon=CortexNL_D1NL
    SynapseMon2=DA_D1NL
    DAind = np.where(SpikeMon2.t > (defaultclock.t - window)) 
    SynapseMon.w = SynapseMon.w * 0.999
    if len(DAind[0]) > 5: # was DA released?
        CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
        PresynapticInd = []
        for i in range(0,len(np.unique(CortexIndex))):
            pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
            PresynapticInd.extend(pre[0])
        PostsynapticInd = np.where(SynapseMon.j[PresynapticInd])  
        NotPostsynapticInd = np.delete(range(0,len(SynapseMon.j)),PostsynapticInd[0],0)
        #PostsynapticNeuron = np.unique(PostsynapticInd) # DA firing is broad, but need this if I want to figure out which neurons recieve both
        #DA ... # DA firing is broad, but need this if I want to figure out which neurons recieve both
        SynapseMon.w[PostsynapticInd[0]] += 3*(SynapseMon.traceCon[PostsynapticInd[0]] * mean(SynapseMon2.traceCon))        
        SynapseMon.w[NotPostsynapticInd[0]] += 3*(SynapseMon.traceConPost[NotPostsynapticInd[0]] * mean(SynapseMon2.traceCon))
        SynapseMon.w = clip(SynapseMon.w, 0, wmax)

@network_operation(dt=window)
def DA_LTP3():
    SpikeMon=CortexWLspikes
    SpikeMon2=DASpikes
    SynapseMon=CortexWL_D1WL
    SynapseMon2=DA_D1WL
    DAind = np.where(SpikeMon2.t > (defaultclock.t - window)) 
    SynapseMon.w = SynapseMon.w * 0.999
    if len(DAind[0]) > 5: # was DA released?
        CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
        PresynapticInd = []
        for i in range(0,len(np.unique(CortexIndex))):
            pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
            PresynapticInd.extend(pre[0])
        PostsynapticInd = np.where(SynapseMon.j[PresynapticInd])  
        NotPostsynapticInd = np.delete(range(0,len(SynapseMon.j)),PostsynapticInd[0],0)
        #PostsynapticNeuron = np.unique(PostsynapticInd) # DA firing is broad, but need this if I want to figure out which neurons recieve both
        #DA ... # DA firing is broad, but need this if I want to figure out which neurons recieve both
        SynapseMon.w[PostsynapticInd[0]] += 3*(SynapseMon.traceCon[PostsynapticInd[0]] * mean(SynapseMon2.traceCon))        
        SynapseMon.w[NotPostsynapticInd[0]] += 3*(SynapseMon.traceConPost[NotPostsynapticInd[0]] * mean(SynapseMon2.traceCon))
        SynapseMon.w = clip(SynapseMon.w, 0, wmax)
#%% Run and analyze
       
def frange(start, stop, step):
 i = start
 while i < stop:
     yield i
     i += step
def calculate_FR(SpikeMon,binSize=100*ms,timeWin=duration):
    allBin = []
    for i in frange(360,int(ceil(timeWin/ms)),binSize/ms):
         binz = len(np.where((SpikeMon.t > i*ms)*(SpikeMon.t < (i*ms+binSize)))[0])/binSize
         allBin.append(binz)
    FR = np.mean(allBin)
    return FR, allBin       
       
sequence = 0
if sequence == 1:  # reproduce figure 3 in humphries et al., 2006      
    ThalamusL.I = 0
    ThalamusNL.I = 0
    ThalamusWL.I = 0
    DA.I = 0
    CortexL.I = 0
    InhInter.I = 5
    run(duration,report='text')
    ThalamusL.I = 0
    ThalamusNL.I = 0
    ThalamusWL.I = 0
    DA.I = 0
    CortexL.I = 50
    CortexNL.I = 0
    InhInter.I = 5
    run(duration,report='text')
    ThalamusL.I = 0
    ThalamusNL.I = 0
    ThalamusWL.I = 0
    DA.I = 0
    CortexL.I = 50
    CortexNL.I = 60
    InhInter.I = 5
    run(duration,report='text')
    
    Lfr,Lbin = calculate_FR(ActionSpikes,binSize=100*ms,timeWin=6)
    NLfr,NLbin = calculate_FR(NoActionSpikes,binSize=100*ms,timeWin=6)
    WLfr,WLbin = calculate_FR(WrongActionSpikes,binSize=100*ms,timeWin=6)
    figure()
    plot(range(0,5700,100),Lbin,'r')
    plot(range(0,5700,100),NLbin,'b')
    plot(range(0,5700,100),WLbin,'g')
    xlabel('Time(ms)')



popFiring = 0
if popFiring == 1: # reproduce figure 2 in Humphries et al., 2006
   #CortexL.I = 15
   #CortexNL.I = 15
   #CortexWL.I = 15   
   run (duration,report='text') 
    
   STN_FR = calculate_FR(STNspikes) # divide by 3 because STN is actually 30 neurons 
   GPe_FR = calculate_FR(GPeSpikes)
   SNr_FR = calculate_FR(SNrLspikes)
   
   print STN_FR[0]/3/2
   print 'STN'
   print GPe_FR[0]/2
   print 'GPe'
   print SNr_FR[0]/2
   print 'SNr'
   
   STNall = np.array([14.618,13.013,14.024,20.206,13.85,14.824,15.69,12.855,16.5189,13.46,13.69])
   STNemp = 10
   SNrAll = np.array([21.319,17.072,20.4,13.57,21.969,25.59,23.7,17.103,24.948,19.40,21.876])
   SNrEmp = 26
   GPeAll = np.array([31.12,26.56,29.154,29.237,29.711,33.577,34.587,26.38144,36.989,27.649,27.98])
   GPeEmp =29

   def mean_confidence_interval(data, confidence=0.95):
        a = 1.0*np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
        return m, m-h, m+h    
    
   STN_95 = mean_confidence_interval(STNall)
   STN_95 = STN_95[2]-STN_95[1]
   SNr_95 = mean_confidence_interval(SNrAll)
   SNr_95 = SNr_95[2]-SNr_95[1]
   GPe_95 = mean_confidence_interval(GPeAll)
   GPe_95 = GPe_95[2]-GPe_95[1]
   
   import matplotlib.pyplot as plt
   figure()
   plt.legend('STN','GPe','SNr')
   plt.plot(1,np.mean(STNall),'ob')
   plt.plot(2,np.mean(GPeAll),'ob')
   plt.plot(3,np.mean(SNrAll),'ob')
   plt.errorbar(1,np.mean(STNall),STN_95,capsize=10, elinewidth=4)
   plt.errorbar(2,np.mean(GPeAll),GPe_95,capsize=10, elinewidth=4)
   plt.errorbar(3,np.mean(SNrAll),SNr_95,capsize=10, elinewidth=4)
   plt.plot(1,STNemp,'g+',mew=10, ms=20)
   plt.plot(2,GPeEmp,'g+', mew=10, ms=20)
   plt.plot(3,SNrEmp,'g+', mew=10, ms=20)
   plt.errorbar(1,STNemp,(14-8)/2,capsize=11, elinewidth=2)
   plt.errorbar(2,GPeEmp,(36-26)/2,capsize=11, elinewidth=2)
   plt.errorbar(3,SNrEmp,(32-21)/2,capsize=11, elinewidth=2)
   xlim(0,4)
   ylim(8,40)
   title('Mean Tonic Firing Rates')

learnAction = 1
if learnAction == 1:
   ThalamusL.I = 0
   ThalamusNL.I = 0
   ThalamusWL.I = 0
   CortexL.I = 2
   CortexNL.I = 2
   CortexWL.I = 2
   run(duration,report='text')
   
   Lfr,Lbin = calculate_FR(ActionSpikes,binSize=100*ms,timeWin=duration/second)
   NLfr,NLbin = calculate_FR(NoActionSpikes,binSize=100*ms,timeWin=duration/second)
   WLfr,WLbin = calculate_FR(WrongActionSpikes,binSize=100*ms,timeWin=duration/second)
   figure()
   plot(range(0,int(duration/ms)-360,100),Lbin,'r')
   plot(range(0,int(duration/ms)-360,100),NLbin,'b')
   plot(range(0,int(duration/ms)-360,100),WLbin,'g')
   xlabel('Time(ms)')
   
   a,SNrLbin = calculate_FR(SNrLspikes,binSize=100*ms,timeWin=duration/second)
   b,SNrNLbin = calculate_FR(SNrNLspikes,binSize=100*ms,timeWin=duration/second)
   c,SNrWLbin = calculate_FR(SNrWLspikes,binSize=100*ms,timeWin=duration/second)
   figure()
   plot(range(0,int(duration/ms)-360,100),SNrLbin,'r')
   plot(range(0,int(duration/ms)-360,100),SNrNLbin,'b')
   plot(range(0,int(duration/ms)-360,100),SNrWLbin,'g')
   xlabel('Time(ms)')
   title('SNr FR')

print np.str(len(ActionSpikes.t)) + '...rewarded action'
print np.str(len(NoActionSpikes.t)) + '...unrewared action'
print np.str(len(WrongActionSpikes.t)) + '...unrewared action'

#%%
def visualise_connectivity(S):
    Ns = len(S.source)
    Nt = len(S.target)
    figure(figsize=(10, 4))
    subplot(121)
    plot(zeros(Ns), arange(Ns), 'ok', ms=10)
    plot(ones(Nt), arange(Nt), 'ok', ms=10)
    for i, j in zip(S.i, S.j):
        plot([0, 1], [i, j], '-k')
    xticks([0, 1], ['Source', 'Target'])
    ylabel('Neuron index')
    xlim(-0.1, 1.1)
    ylim(-1, max(Ns, Nt))
    subplot(122)
    plot(S.i, S.j, 'ok')
    xlim(-1, Ns)
    ylim(-1, Nt)
    xlabel('Source neuron index')
    ylabel('Target neuron index')
    
#visualise_connectivity(NoPress_Press)    
#visualise_connectivity(SNr_Action)

figure()
plot(ActionSpikes.t/ms,ActionSpikes.i,'xb')
plot(NoActionSpikes.t/ms,NoActionSpikes.i,'or')
plot(WrongActionSpikes.t/ms,WrongActionSpikes.i,'.g')
plot(CortexLSpikes.t/ms,np.ones([1,len(CortexLSpikes.t)])[0],'c*')
ylim(-1,2)
xlabel('Time(ms)')
ylabel('Neuron Number')
title('Brainstem Spikes')
