# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 14:31:28 2016

@author: Elijah
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 09:42:16 2016

@author: Elijah
"""

#%%
# TO DO

# make two seperate regions of the thalamus. One will decide an action the other will feed to the cortex. 
#%%


from brian2 import *
#prefs.codegen.target = 'weave'
start_scope()
duration = 1000*ms#1000000*ms

recordz = 1
plotz = 1


#%%
# Parameters

n = 50 # number of neurons
w2 = n # for connectivity

# thalamo-cortical Iz
c = -66#-66*mV
b= 0.2#0.2/ms 
a = 0.02#0.02/ms
d = 8#8*mV/ms

# cortical
taum = 5*ms#2*ms#5*ms#10*ms
taupsp = 0.325*ms#0.325*ms
weight = 4.86*mV

synfire = 1
recurrent = 0

# DA 
DAltpRate = 100 # the factor synapses that got DA and cortical input will be potentiated by
DAltdRate = 0.6 # factor to decrease D2 activity by
ltpDecay = 0.995 # amount LTP decreases with each cortical spike
window = 60*ms # amount of time for DA integration.. from Yagishita et al.,2014; 0.3-2s

# sensory Threshold used for knowing if tone is on
sensoryThreshold = -30 # taking mean activity at the moment hence the negative

#%% Equations 

# Missing the (1/.taum)
#eqs =  '''
#dv/dt = (((0.04)*v**2+(5)*v+140-u+I)+x)/ms : 1
#du/dt = a*(b*v-u)/ms : 1
#dx/dt = ((-x+y/mV)*(1./taupsp)) : 1
#dy/dt = -y*(1./taupsp)+25.27*mV/ms+(39.24*mV/ms**0.5)*xi : volt
#I : 1
#'''

# equations similar to original brian synfire example... 
# The groups are limited to about 4 and span larger distances 
eqs =  '''
dv/dt = ((((0.04)*v**2+(5)*v+140-u+I)+x)*((1./taum)*ms))/ms : 1
du/dt = a*(b*v-u)/ms : 1
dx/dt = ((-x+y/mV)*(1./taupsp)) : 1
dy/dt = -y*(1./taupsp)+25.27*mV/ms+(39.24*mV/ms**0.5)*xi : volt
I : 1
'''

eqs2 = '''
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = a*(b*v-u)/ms : 1
I : 1
'''

reset='''
v=c
u += d
'''

MSNeqs = '''
dv/dt = (((v+80)*(v+25)-u+I)/ms)/50 : 1
du/dt = 0.01*(-20*(v+150))/ms : 1
I : 1
'''

MSNreset = '''
v = -55
u+=150
'''

#%% Neuron groups
n_groups = 10
group_size = 80 # This will result in 800 neurons which is what Laje and Buonomano 2013 used in their RNN
w3 = n_groups*group_size

# Poisson Input
P_STN = PoissonGroup(n, np.arange(n)*Hz + 10*Hz) # input to STN neurons

# Actions
LeverPress = NeuronGroup(10,eqs2,threshold='v>60',reset=reset,method='euler')
WrongLeverPress = NeuronGroup(10,eqs2,threshold='v>60',reset=reset,method='euler')
NoLeverPress = NeuronGroup(10,eqs2,threshold='v>20',reset=reset,method='euler')

Cortex = NeuronGroup(N=n_groups*group_size, model=eqs,
threshold='v>30', reset=reset, refractory=2*ms,
method='euler')

# DA neurons
DA = NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')

DA.v = c
DA.u = b*c

# medium spiny neurons 
D1=NeuronGroup(n,MSNeqs,threshold='v > 40',
                    reset=MSNreset,
                    method='euler')
D1.v = c
D1.u = b*c

D2=NeuronGroup(n,MSNeqs,threshold='v > 40',
                    reset=MSNreset,
                    method='euler')
D2.v = c
D2.u = b*c

#GPe neurons
GPe=NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')

GPe.v = c
GPe.u = b*c

# Hypothalamus
Hypothalamus=NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')

Hypothalamus.v = c
Hypothalamus.u = b*c

# SNr neurons
SNr=NeuronGroup(n,eqs2,threshold='v>30', reset=reset,method='euler')

SNr.v = c
SNr.u = b*c

# STN neurons
STN=NeuronGroup(n,eqs2,threshold='v>20',reset=reset,method='euler') #ampaEqs,threshold='v>20*mV',reset=reset4,method='euler')

STN.v = c
STN.u = b*c

# Thalamic neurons
Thalamus =NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')

Thalamus.v = c
Thalamus.u = b*c

#%% The network structure
# Synapses

# STDP variables
taupre = taupost = 20*ms
wmax = 1 # weights can go from 0 to 1
Apre = 0.01
Apost = -Apre*taupre/taupost*1.05
wScale = 1 # amount to scale weights by

# MSN plasticity
traceTau = traceTauPost = 1200*ms  # this will change integration window
traceConstant = 0.01 # this will change maximum weight change
tracePostConstant = -traceConstant*traceTau/traceTauPost*1.05

# STDP eqs
stdp_eqs = '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             '''
onpre='''
             v_post += w*wScale
             apre += Apre
             w = clip(w+apost, 0, wmax)
             '''
onpreNeg='''
             v_post -= w*wScale
             apre += Apre
             w = clip(w+apost, 0, wmax)
             '''   
onpost='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             '''
             
             
DAstdp = '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             dtraceCon/dt = -traceCon/traceTau : 1 (event-driven)
             dtraceConPost/dt = -traceConPost/traceTauPost : 1 (event-driven)
             '''   
onpreDA_D1='''
             v_post *= wmax+w*wScale
             apre += Apre
             w = clip(w+apost, 0, wmax)
             traceCon += traceConstant
             '''
onpreDA_D2='''
             v_post *= wmax-(w*wScale)
             w = clip(w+apost, 0, wmax)
             traceCon += traceConstant
             apre += Apre
             '''                                                  
onpostDA='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             traceConPost += tracePostConstant 
         '''

        
MSNstdp=  '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             dtraceCon/dt = -traceCon/traceTau : 1 (event-driven)
             dtraceConPost/dt = -traceConPost/traceTauPost : 1 (event-driven)
             '''    
MSNpre='''
             v_post += w
             apre += Apre
             w = clip(w+apost, 0, wmax)
             traceCon += traceConstant
             '''

MSNpost ='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             traceConPost += tracePostConstant 
             '''

################################################################
# Poisson Input
#P_STN_STN = Synapses(P_STN,STN, on_pre = 'v+=10')
#P_STN_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) 
#P_STN_STN.delay = 10*ms

################################################################
# Excitatory Connections
Cortex_D1 = Synapses(Cortex, D1, MSNstdp,on_pre=MSNpre,on_post=MSNpost)
Cortex_D1.connect(j='k for k in range(i-w3, i+w3) if rand()<0.5', skip_if_invalid=True)
#Cortex_D1.connect(condition='i!=j',p=0.5,skip_if_invalid=True)
Cortex_D1.delay = 10*ms
Cortex_D1start = rand(len(Cortex_D1.i))
Cortex_D1.w = Cortex_D1start

Cortex_D2 = Synapses(Cortex, D2, MSNstdp,on_pre=MSNpre,on_post=MSNpost)
Cortex_D2.connect(j='k for k in range(i-w3, i+w3) if rand()<0.5', skip_if_invalid=True)
#Cortex_D2.connect(condition='i!=j',p=0.5,skip_if_invalid=True)
Cortex_D2.delay = 10*ms
Cortex_D2.w = rand(len(Cortex_D2.i))

Cortex_DA = Synapses(Cortex,DA,stdp_eqs,on_pre=onpre,on_post=onpost)  
Cortex_DA.connect(j='k for k in range(i-w3, i+w3) if rand()<0.2', skip_if_invalid=True) # don't think this should be too high or get strong motivation for action all the time 
#Cortex_DA.connect(condition='i!=j',p=0.6)
Cortex_DA.delay = 10*ms
Cortex_DA.w = rand(len(Cortex_DA.i))

Cortex_LeverPress = Synapses(Cortex,LeverPress,stdp_eqs,on_pre=onpre,on_post=onpost)    
Cortex_LeverPress.connect(j='k for k in range(i-w3, i+w3) if rand()<0.2', skip_if_invalid=True)  
#Cortex_LeverPress.connect(condition='i!=j',p=0.5,skip_if_invalid=True)
Cortex_LeverPress.delay = 10*ms
Cortex_LeverPress.w = rand(len(Cortex_LeverPress.i))

Cortex_WrongLeverPress = Synapses(Cortex,WrongLeverPress,stdp_eqs,on_pre=onpre,on_post=onpost)    
Cortex_WrongLeverPress.connect(j='k for k in range(i-w3, i+w3) if rand()<0.2', skip_if_invalid=True)  
#Cortex_LeverPress.connect(condition='i!=j',p=0.5,skip_if_invalid=True)
Cortex_WrongLeverPress.delay = 10*ms
Cortex_WrongLeverPress.w = rand(len(Cortex_WrongLeverPress.i))

Cortex_NoLeverPress = Synapses(Cortex,NoLeverPress,stdp_eqs,on_pre=onpre,on_post=onpost) 
Cortex_NoLeverPress.connect(j='k for k in range(i-w3, i+w3) if rand()<0.3', skip_if_invalid=True)     
#Cortex_NoLeverPress.connect(condition='i!=j',p=0.5,skip_if_invalid=True)
Cortex_NoLeverPress.delay = 10*ms
Cortex_NoLeverPress.w = rand(len(Cortex_NoLeverPress.i))

# Hyperdirect pathway
Cortex_STN = Synapses(Cortex,STN,stdp_eqs,on_pre=onpre,on_post=onpost) 
Cortex_STN.connect(j='k for k in range(i-w3, i+w3) if rand()<0.5', skip_if_invalid=True)  
Cortex_STN.delay = 5*ms # tis is in order to make the hyperdirect pathway quicker than the direct pathway
Cortex_STN.w = rand(len(Cortex_STN.i))

# Excitatory Hypothalamic projections to DA
S_Hypo_DA = Synapses(Hypothalamus,DA,stdp_eqs,on_pre=onpre,on_post=onpost) 
S_Hypo_DA.connect(j='k for k in range(i-w2, i+w2) if rand()<0.7', skip_if_invalid=True)  
#S_Hypo_DA.connect(condition='i!=j',p=0.8,skip_if_invalid=True)
S_Hypo_DA.delay = 5*ms
S_Hypo_DA.w = rand(len(S_Hypo_DA.i))

# Hyperdirect pathway
STN_SNr = Synapses(STN, SNr,stdp_eqs,on_pre=onpre,on_post=onpost)
STN_SNr.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)  
STN_SNr.delay = 2*ms # This is just matching the STN to GPe
STN_SNr.w = rand(len(STN_SNr.i))

# Excitatory STN to GPe
S_STN_GPe = Synapses(STN,GPe,stdp_eqs,on_pre=onpre,on_post=onpost)
S_STN_GPe.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True)  
#S_STN_GPe.connect(condition='i!=j',p=0.2,skip_if_invalid=True)
S_STN_GPe.delay = 2*ms
S_STN_GPe.w = rand(len(S_STN_GPe.i))

#Thalamic input to cortex; excitatory 
Thal_Cortex = Synapses(Thalamus,Cortex,stdp_eqs,on_pre=onpre,on_post=onpost)
Thal_Cortex.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True)  
Thal_Cortex.delay = 2*ms
Thal_Cortex.w = rand(len(Thal_Cortex.i))

# Thalamus to Striatum; excitatory
Thal_D1 = Synapses(Thalamus,D1,stdp_eqs,on_pre=onpre,on_post=onpost)
Thal_D1.connect(j='k for k in range(i-w2, i+w2) if rand()<0.3', skip_if_invalid=True)  
#Thal_D1.connect(condition='i!=j',p=0.3,skip_if_invalid=True) # used lower probability than cortex
Thal_D1.delay = 10*ms # used same delay as cortex
Thal_D1.w = rand(len(Thal_D1.i))

Thal_D2 = Synapses(Thalamus,D2,stdp_eqs,on_pre=onpre,on_post=onpost)
Thal_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.3', skip_if_invalid=True)  
#Thal_D2.connect(condition='i!=j',p=0.3,skip_if_invalid=True) # used lower probability than cortex
Thal_D2.delay = 10*ms # used same delay as cortex
Thal_D2.w = rand(len(Thal_D2.i))

################################################################
# Inhibitory connections 

# BG input to GPe; inhibitory
BG_GPe = Synapses(D2,GPe,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
BG_GPe.connect(j='k for k in range(i-w2, i+w2) if rand()<0.34', skip_if_invalid=True)  
#BG_GPe.connect(condition='i!=j',p=0.34,skip_if_invalid=True)
BG_GPe.delay = 5*ms
BG_GPe.w = rand(len(BG_GPe.i))

# BG input to SNr; inhibitory
BG_SNr = Synapses(D1,SNr,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
BG_SNr.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#BG_SNr.connect(condition='i!=j',p=1,skip_if_invalid=True)
BG_SNr.delay = 4*ms
BG_SNr.w = rand(len(BG_SNr.i))

# GPe to STN; Inhibitory
GPe_STN = Synapses(GPe,STN,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
GPe_STN.connect(j='k for k in range(i-w2, i+w2) if rand()<0.4', skip_if_invalid=True) 
#GPe_STN.connect(condition='i!=j',p=0.4,skip_if_invalid=True)
GPe_STN.delay = 4*ms
GPe_STN.w = rand(len(GPe_STN.i))

# SNr to thalamic; Inhibitory
SNr_Thal = Synapses(SNr,Thalamus,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
SNr_Thal.connect(j='k for k in range(i-w2, i+w2) if rand()<0.3', skip_if_invalid=True) 
#SNr_Thal.connect(condition='i!=j',p=0.3,skip_if_invalid=True)
SNr_Thal.delay = 4*ms
SNr_Thal.w = rand(len(SNr_Thal.i))

SNr_LeverPress = Synapses(SNr,LeverPress,stdp_eqs,on_pre=onpreNeg,on_post=onpost) # This will be similar to connectivity to PPN; Humphries, M. D., Gurney, K., & Prescott, T. J. (2007). Is there a brainstem substrate for action selection?. Philosophical Transactions of the Royal Society of London B: Biological Sciences, 362(1485), 1627-1639.
SNr_LeverPress.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
#SNr_Thal.connect(condition='i!=j',p=0.3,skip_if_invalid=True)
SNr_LeverPress.delay = 5*ms
SNr_LeverPress.w = rand(len(SNr_LeverPress.i))

SNr_WrongLeverPress = Synapses(SNr,WrongLeverPress,stdp_eqs,on_pre=onpreNeg,on_post=onpost) # This will be similar to connectivity to PPN; Humphries, M. D., Gurney, K., & Prescott, T. J. (2007). Is there a brainstem substrate for action selection?. Philosophical Transactions of the Royal Society of London B: Biological Sciences, 362(1485), 1627-1639.
SNr_WrongLeverPress.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
#SNr_Thal.connect(condition='i!=j',p=0.3,skip_if_invalid=True)
SNr_WrongLeverPress.delay = 5*ms
SNr_WrongLeverPress.w = rand(len(SNr_WrongLeverPress.i))

SNr_NoLeverPress = Synapses(SNr,NoLeverPress,stdp_eqs,on_pre=onpreNeg,on_post=onpost) # This will be similar to connectivity to PPN; Humphries, M. D., Gurney, K., & Prescott, T. J. (2007). Is there a brainstem substrate for action selection?. Philosophical Transactions of the Royal Society of London B: Biological Sciences, 362(1485), 1627-1639.
SNr_NoLeverPress.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
#SNr_Thal.connect(condition='i!=j',p=0.3,skip_if_invalid=True)
SNr_NoLeverPress.delay = 5*ms
SNr_NoLeverPress.w = rand(len(SNr_NoLeverPress.i))

################################################################
# Collaterals and recurrents

# Action to Action
# this is building a winner take all circuit
LeverPress_NoPress = Synapses(LeverPress,NoLeverPress,on_pre='v-=5') # These are inhibitory with longer delays
LeverPress_NoPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_NoPress.delay = 5*ms

LeverPress_NoPressEx = Synapses(LeverPress,NoLeverPress,on_pre='v+=5') # excitatory with lower delays 
LeverPress_NoPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_NoPressEx.delay = 2*ms

LeverPress_WrongPress = Synapses(LeverPress,WrongLeverPress,on_pre='v-=5')
LeverPress_WrongPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_WrongPress.delay = 5*ms

LeverPress_WrongPressEx = Synapses(LeverPress,WrongLeverPress,on_pre='v+=5')
LeverPress_WrongPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_WrongPressEx.delay = 2*ms

# this is building a winner take all circuit
WrongPress_NoPress = Synapses(WrongLeverPress,NoLeverPress,on_pre='v-=5')
WrongPress_NoPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_NoPress.delay = 5*ms

WrongPress_NoPressEx = Synapses(WrongLeverPress,NoLeverPress,on_pre='v+=5')
WrongPress_NoPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_NoPressEx.delay = 2*ms

WrongPress_Press = Synapses(WrongLeverPress,LeverPress,on_pre='v-=5')
WrongPress_Press.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_Press.delay = 5*ms

WrongPress_PressEx = Synapses(WrongLeverPress,LeverPress,on_pre='v+=5')
WrongPress_PressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_PressEx.delay = 2*ms


NoPress_Press = Synapses(NoLeverPress,LeverPress,on_pre='v-=10')
NoPress_Press.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#NoPress_Press.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_Press.delay = 5*ms
#NoPress_Press.w = rand(len(NoPress_Press.i))

NoPress_PressEx = Synapses(NoLeverPress,LeverPress,on_pre='v+=5')
NoPress_PressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#NoPress_Press.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_PressEx.delay = 2*ms

NoPress_WrongPress = Synapses(NoLeverPress,WrongLeverPress,on_pre='v-=10')
NoPress_WrongPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#NoPress_Press.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_WrongPress.delay = 5*ms
#NoPress_Press.w = rand(len(NoPress_Press.i))

NoPress_WrongPressEx = Synapses(NoLeverPress,WrongLeverPress,on_pre='v+=5')
NoPress_WrongPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#NoPress_Press.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_WrongPressEx.delay = 5*ms

# Inhibitory SNr collaterals 
S_SNr_SNr = Synapses(SNr,SNr,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
S_SNr_SNr.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True) 
#S_SNr_SNr.connect(condition='i!=j',p=0.5,skip_if_invalid=True)
S_SNr_SNr.delay = 1*ms
S_SNr_SNr.w = rand(len(S_SNr_SNr.i))

# Inhibitory Striatal collaterals
D1_D1 = Synapses(D1,D1,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
D1_D1.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
#D1_D1.connect(condition='i!=j',p=0.1,skip_if_invalid=True)
D1_D1.delay = 1*ms
D1_D1.w = rand(len(D1_D1.i))

D2_D2 = Synapses(D2,D2,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
D2_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.1', skip_if_invalid=True) 
#D2_D2.connect(condition='i!=j',p=0.1,skip_if_invalid=True)
D2_D2.delay = 1*ms
D2_D2.w = rand(len(D2_D2.i))

# synfire connections
if synfire == 1:
    Cortex_Cortex = Synapses(Cortex, Cortex, on_pre='y+=weight')
    Cortex_Cortex.connect(j='k for k in range((int(i/group_size)+1)*group_size, (int(i/group_size)+2)*group_size) '
    'if i<N_pre-group_size')
    #Cortex_Cortex.delay = 10*ms
if recurrent == 1:
    Cortex_Cortex = Synapses(Cortex,Cortex,on_pre='y+=weight') # this keeps similar amounts of acitity in the network
    Cortex_Cortex.connect(condition='i!=j',p=0.1,skip_if_invalid=True) # probability used in Laje & Buonomano, 2013
    Cortex_Cortex.delay = 1*ms
    Cortex_Cortex.w = rand(len(Cortex_Cortex.i))
    
##################################################################
# Modulatory 

# DA projections to D1
DA_D1 = Synapses(DA,D1,DAstdp,on_pre=onpreDA_D1,on_post=onpostDA)
DA_D1.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
#DA_D1.connect(condition='i!=j',p=0.75,skip_if_invalid=True)  #p=0.1) # I think I'd rather have a broad range of DA with synapse specific plasticity at places with glut and DA
DA_D1.delay = 5*ms
DA_D1.w = rand(len(DA_D1.i))

# DA projections to D2
DA_D2 = Synapses(DA,D2,DAstdp,on_pre=onpreDA_D2,on_post=onpostDA)
DA_D2.connect(j='k for k in range(i-w2, i+w2) if rand()<0.75', skip_if_invalid=True) 
#DA_D2.connect(condition='i!=j',p=0.75,skip_if_invalid=True) # p0.1
DA_D2.delay = 5*ms
DA_D2.w = rand(len(DA_D2))

#%% Network Operators

 #######3
# The problem is I am currently using the indicies of the presynaptic neurons, Which are not 
# the same for the two populations of neurons. Need to find out when postsynaptic neurons got presynaptic input


#@network_operation(dt=10*ms)
#def DA_LTD():
#    #CortexIndex = Cortex_D2.i[Cortex_D2.t > defaultclock.t-window]    
#    CortexIndex = cortexSpikes.i[cortexSpikes.t > defaultclock.t-window]  # index for the active neuron
#    DAindex = DAspikes.i[DAspikes.t > defaultclock.t-window]
#    bothIndex = intersect1d(array(CortexIndex),array(DAindex))
#    synapseIndexArray = Cortex_D2.i[0:np.size(Cortex_D2.i)]
#    synapseIndexMat = []
#    for j in range(0,np.size(bothIndex)):
#        if np.size(bothIndex) > 0:
#            synapseIndex = np.where(synapseIndexArray==bothIndex[j])[0]
#            synapseIndexMat.extend(synapseIndex)    
#    #print synapseIndexMat
#    #Cortex_D2.w_[synapseIndexMat] *= DAltdRate
#    Cortex_D2.w_[synapseIndexMat] = clip(Cortex_D2.w_[synapseIndexMat]+(-Cortex_D2.apost_[synapseIndexMat]*DAltdRate),0,wmax)
  
#@network_operation(dt=10*ms)
## the problem with this is that I am not looking at presynaptic and post synaptic neurons
##consider using Cortex_D1.source and DA_D1.source
#def DA_LTP():
#    CortexIndex = cortexSpikes.i[cortexSpikes.t > defaultclock.t-window]  # index for the active neuron
#    DAindex = DAspikes.i[DAspikes.t > defaultclock.t-window]
#    bothIndex = intersect1d(array(CortexIndex),array(DAindex))
#    synapseIndexArray = Cortex_D1.i[0:np.size(Cortex_D1.i)]
#    synapseIndexMat = []
#    for j in range(0,np.size(bothIndex)):
#        if np.size(bothIndex) > 0:
#            synapseIndex = np.where(synapseIndexArray==bothIndex[j])[0]
#            synapseIndexMat.extend(synapseIndex)    
#    #print synapseIndexMat
#    #Cortex_D1.w_[synapseIndexMat] *= DAltpRate #1/(1+math.exp(DAltpRate))                  
#    Cortex_D1.w_[synapseIndexMat] = clip(Cortex_D1.w_[synapseIndexMat]+(Cortex_D1.apost_[synapseIndexMat]*DAltpRate),0,wmax)       
#    #Cortex_D1.w_[synapseIndexMat] = clip(Cortex_D1.w_[synapseIndexMat]*DAltpRate),0,wmax)       

@network_operation(dt=10*ms)
def DA_LTD():
    indz = np.where(DA_D1.t > (defaultclock.t - 10*ms)) #Cortex_D2_DAtrace.t > defaultclock.t-10*ms)
    if len(indz) > 2:
       Cortex_D2.w =- Cortex_D2.traceCon * DA_D2.traceCon
       #Cortex_D2_DAtrace = 0 # resetting variables to save memroy

@network_operation(dt=10*ms)
def DA_LTP():
    indz = np.where(DA_D2.t > (defaultclock.t - 10*ms)) #Cortex_D1_DAtrace.t > defaultclock.t-10*ms)
    if len(indz) > 2:
       Cortex_D1.w =+ Cortex_D1.traceCon * DA_D1.traceCon
       #Cortex_D1_DAtrace  # resetting variables to save memory
    
@network_operation(dt=10*ms) # here the reward rule will be give reward if lever press occurs druing tone
def DA_Reward():
    if mean(Thalamus.v) > sensoryThreshold:
       actionSelection = np.array([len(ActionSpikes.t[ActionSpikes.t > defaultclock.t-10*ms]),len(NoActionSpikes.t[NoActionSpikes.t > defaultclock.t-10*ms]),len(WrongActionSpikes.t[WrongActionSpikes.t > defaultclock.t-10*ms])])      
       if np.where(actionSelection == np.max(actionSelection))<1:
           Hypothalamus.I   -= -0.1 #0.01*(len(ActionSpikes.t[ActionSpikes.t > defaultclock.t-10*ms])) #decrease motivation
           DA.v += 2 #len(ActionSpikes.t[ActionSpikes.t > defaultclock.t-10*ms]) # increase reward

@network_operation(dt=500*ms) # update the sensory input every 100ms
def sensInput():
    if Thalamus.I[0] > 0:
       Thalamus.I[:] = 0
    else: 
        Thalamus.I[:] = 50
    #Thalamus.I = 100*(sin(defaultclock.t/ms)+0.5)     



#@network_operation(dt=100*ms)
#def recordStartWeights():
#    if defaultclock.t < 120*ms:
#
#       Cortex_D1start.w = Cortex_D1.w # record starting weights
#       Cortex_D1startMon = StateMonitor(Cortex_D1start,'w',record=True)
#       print 'here'

    
#%% Record the spikes

ActionSpikes = SpikeMonitor(LeverPress)
NoActionSpikes = SpikeMonitor(NoLeverPress)
WrongActionSpikes = SpikeMonitor(WrongLeverPress)

# population activity
ActionPop = PopulationRateMonitor(LeverPress)
NoActionPop = PopulationRateMonitor(NoLeverPress)
WrongActionPop = PopulationRateMonitor(WrongLeverPress)
D1pop = PopulationRateMonitor(D1)
HypoPop  = PopulationRateMonitor(Hypothalamus)
ThalamusPop = PopulationRateMonitor(Thalamus)
DApop = PopulationRateMonitor(DA)

# Weights
#D1weights = StateMonitor(Cortex_D1,'w',record=range(2))
CortexD1change = Cortex_D1start-Cortex_D1.w

from matplotlib.pyplot import pcolormesh
#pcolormesh(1,range(0,len(CortexD1change)),CortexD1change)
#import numpy as np
#np.repeat(CortexD1change,1)


if recordz == 1:    
    ActionTrace = StateMonitor(LeverPress,('v'),record=True)  
    NoActionTrace = StateMonitor(NoLeverPress,('v'),record=True)
    WrongActionTrace = StateMonitor(WrongLeverPress,('v'),record=True)
    
    cortexSpikes = SpikeMonitor(Cortex)
    cortexTrace = StateMonitor(Cortex, 'v', record=True)
    
    D1_Spikes = SpikeMonitor(D1)
    D1_Trace = StateMonitor(D1, ('v'), record=True)
    
    D2_Spikes = SpikeMonitor(D2)
    D2_Trace = StateMonitor(D2, ('v'), record=True)
    
    DAspikes = SpikeMonitor(DA)
    DAtrace = StateMonitor(DA,'v',record=True)
    
    HypoSpikes  = SpikeMonitor(Hypothalamus)
    HypoTrace = StateMonitor(Hypothalamus, ('v'), record=True)
    
    hypoMon = StateMonitor(Hypothalamus,('I'),record=True)   
    
    SNrSpikes = SpikeMonitor(SNr)
    SNrTrace = StateMonitor(SNr, ('v'), record=True)
    
    thalamusSpikes = SpikeMonitor(Thalamus)
    thalamusTrace = StateMonitor(Thalamus, ('v'), record=True)
    
    GPeSpikes = SpikeMonitor(GPe)
    GPeTrace = StateMonitor(GPe, ('v'), record=True)
    
    STNSpikes = SpikeMonitor(STN)
    STNTrace = StateMonitor(STN, ('v'), record=True)
    
    sensMon = StateMonitor(Thalamus,'I',record=True)
    
    Cortex_D1mon = StateMonitor(Cortex_D1,'w', record=Cortex_D1['i!=j'])
    Cortex_D2mon = StateMonitor(Cortex_D2,'w', record=Cortex_D2['i!=j'])  
    
    # Monitors for Glutamate-DA plasticity
    Cortex_D1sigTrace = StateMonitor(Cortex_D1,'traceCon', record=Cortex_D1['i!=j'])
    Cortex_D1_DAtrace = StateMonitor(DA_D1,'traceCon', record=DA_D1['i!=j'])
        
    Cortex_D2sigTrace = StateMonitor(Cortex_D2,'traceCon', record=Cortex_D2['i!=j'])
    Cortex_D2_DAtrace = StateMonitor(DA_D2,'traceCon', record=DA_D2['i!=j'])

  
    

    

    
#%%
# run https://mail.google.com/mail/u/0/?ui=2&ik=3ce1d7c4a6&view=att&th=158d10a24fb9809e&attid=0.1&disp=safe&realattid=f_iwcms82h1&zw
Hypothalamus.I = 80 # motivation
Thalamus.I = 0 #sensory stimulus
run(duration, report='text')
#Thalamus.I = 100*(sin(defaultclock.t/ms)+0.5) # sensory stimulus
#run(duration, report='text')
#Thalamus.I = 100*(sin(defaultclock.t/ms)+0.5) # sensory stimulus
#run(duration, report='text')
#Thalamus.I = 100*(sin(defaultclock.t/ms)+0.5) # sensory stimulus
#run(duration, report='text')
#Thalamus.I = 100*(sin(defaultclock.t/ms)+0.5) # sensory stimulus
#run(duration, report='text')

#%% Plot
figure()
plot(ActionPop.t/ms,ActionPop.rate/n/Hz)
title('Lever Press Activity')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

figure()
plot(D1pop.t/ms,D1pop.rate/n/Hz)
title('D1 Population Activity')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

figure()
plot(HypoPop.t/ms,HypoPop.rate/n/Hz)
title('Hypo Population Activity')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

figure()
plot(ThalamusPop.t/ms,ThalamusPop.rate/n/Hz)
title('Thalamus Population Activity')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

figure()
plot(DApop.t/ms, DApop.rate/n/Hz)
title('DA Poplation Activity')
xlabel('Time (ms)')
ylabel('Rate (Hz)')



if recordz == 1:
    figure()
    plot(cortexSpikes.t/ms, 1.0*cortexSpikes.i/group_size, '.')
    plot([0, duration/ms], np.arange(n_groups).repeat(2).reshape(-1, 2).T, 'k-')
    ylabel('group number')
    yticks(np.arange(n_groups))
    xlabel('time (ms)')
    show()
    
    figure()
    plot(cortexTrace.t/ms,cortexTrace.v[1])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('Cortical Voltage')
    
    figure()
    plot(D1_Trace.t/ms,D1_Trace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('D1_MSN Voltage')
    
    figure()
    plot(D2_Trace.t/ms,D2_Trace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('D2_MSN Voltage')
    
    figure()
    plot(SNrTrace.t/ms,SNrTrace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('SNr Voltage')
    
    figure()
    plot(thalamusTrace.t/ms,thalamusTrace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('Thalamic Voltage')
    
    figure()
    plot(DAtrace.t/ms,DAtrace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('DA Voltage')
    
    figure()
    plot(HypoTrace.t/ms,HypoTrace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('Hypothalamus Voltage')
    
    figure()
    plot(GPeTrace.t/ms,GPeTrace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('GPe Voltage')
    
    figure()
    plot(STNTrace.t/ms,STNTrace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('STN Voltage')
    
    figure()
    plot(Cortex_D1mon.t/ms, Cortex_D1mon.w[0])
    xlabel('Time(ms)')
    ylabel('Weights')
    title('Cortical D1 type weight')
    
    figure()
    plot(Cortex_D2mon.t/ms, Cortex_D2mon.w[0])
    xlabel('Time(ms)')
    ylabel('Weights')
    title('Cortical D2 type weight')
    
    figure()
    plot(ActionTrace.t/ms,ActionTrace.v[0])
    plot(NoActionTrace.t/ms,NoActionTrace.v[0])
    plot(WrongActionTrace.t/ms,WrongActionTrace.v[0])
    xlabel('Time(ms)')
    ylabel('Voltage(mV)')
    title('Brainstem Voltage')
    
    figure()
    plot(hypoMon.t/ms,hypoMon.I[0])
    xlabel('Time(ms)')
    ylabel('Input Current (amp)')
    title('Hypothalamus "motivation" over time')
    
    checkSTDP = 1
    if checkSTDP == 1:
        figure()
        delta_t = linspace(-50, 50, 100)*ms
        W = where(delta_t<0, Apre*exp(delta_t/taupre), Apost*exp(-delta_t/taupost))
        plot(delta_t/ms, W)
        xlabel(r'$\Delta t$ (ms)')
        ylabel('W')
        ylim(-Apost, Apost)
        axhline(0, ls='-', c='k')
        
        figure()
        delta_t = linspace(-5000, 5000, 100)*ms
        W2 = where(delta_t<0, traceConstant*exp(delta_t/traceTau), tracePostConstant*exp(-delta_t/traceTauPost))
        plot(delta_t/ms, W2)
        xlabel(r'$\Delta t$ (ms)')
        ylabel('W')
        ylim(-tracePostConstant, tracePostConstant)
        axhline(0, ls='-', c='k')
        
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
    
    visualise_connectivity(Cortex_D1)    
    
        
        