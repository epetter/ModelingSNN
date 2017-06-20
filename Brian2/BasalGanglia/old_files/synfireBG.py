# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:54:54 2016

@author: Elijah
"""


#make this model select an action. The action will have some outcome based upon state
# Send outcome to hypothalamus. If outcome decreases activity (e.g. food if hungry), decrease hypothalamus activity
# hypothalamus activity drives actions, actions of the appropriate type inhibit hypothalamic activity
# keep making actions until hypothalamic acitvity is reduced. 

# Could also have sensory inut (VIA THALAMUS) that signal time of reward, or stimuli in environment
# this may allow for development of timing
# motivation will vary in time-dependent manner

# simulate some motivational state with hypothalamic activity. Brainstem actions will inhbit hypothalamus activity
# also feed sensory information into the thalamus. Can just do this through time varying current inputs (e.g. start of tone)

# as of right now all widths are the same as I redefine them 

# Iz neuron network
from brian2 import *
start_scope()

#%%
####################################################################
############ variables #############################################
####################################################################

# Parameters
DAltpRate = 100 # the factor synapses that got DA and cortical input will be potentiated by
DAltdRate = 0.6 # factor to decrease D2 activity by
ltpDecay = 0.995 # amount LTP decreases with each cortical spike
window = 60*ms # amount of time for DA integration.. from Yagishita et al.,2014; 0.3-2s
feedbackConstant = 200 # constant for determining reduction of error signal

# Synfire chains 
taum = 10*ms
taupsp = 0.325*ms
weight = 0*mV#4.86*mV


# IZ neurons
n = 100;
duration = 100*ms
taus = 1*ms
neuronType1 ='tc'
neuronType2 = 'rs' # will change to msn when I get a chance to update equations
neuronType3 = 'rs' # this will be for most regular neurons 

# thalamo-cortical Iz
c = -66*mV#
b= 0.2/ms
a = 0.02/ms
d = 8*mV/ms
# dimensionsless t-c Iz
#c = -66
#b= 0.2
#a = 0.02
#d = 8


# medium spiny neurons... need to adjust equations as well    
msn_c = -66*mV#-55*mV
msn_b = 0.2/ms#-20/ms#-20/ms
msn_a = 0.02/ms#0.01/ms
msn_d = 8*mV/ms#150*mV/ms#150*mV/ms
#https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=128818#tabs-2

#cortical neurons
c3 = -66*mV
b3= 0.2/ms
a3 = 0.1/ms
d3 = 2*mV/ms

# RS neurons for GPe and STN
c4 = -66*mV
b4= 0.2/ms
a4 = 0.02/ms
d4 = 8*mV/ms


a5 = -66
b5 = 0.2
c5 = 0.02
d5 = 8


#%%
####################################################################
############ eqs ###################################################
####################################################################


eqs =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = msn_a*(msn_b*v-u) : volt/second
I : amp
x : metre
y : metre
area : 1
previous_v : volt
'''

eqs2 =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = msn_a*(msn_b*v-u) : volt/second
I : amp
x : metre
'''

#eqs2 = '''
#dv/dt = (1/ms/mV)*v**2+(20/ms)*v+80*mV/ms-u+I/nF : volt
#du/dt = msn_a*(msn_b*(v+80*mV)-u) : volt/second
#I : amp
#x : metre
#'''

#dv/dt = (v+80*mV)*(v+25*mV)/mV/ms-u+I/nF : volt


eqs3 = '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a3*(b3*v-u) : volt/second
I : amp
x : metre
'''

eqs4 = '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a4*(b4*v-u) : volt/second
I : amp
x : metre
'''

#ampaEqs = '''
#dv/dt = (0.04)*v**2+(5)*v+140*mV-u/ms+I+g_ampa/nS : 1
#du/dt = a5*(b5*v-u) : Hz
#I : 1
#x : metre
#dg_ampa/dt = -g_ampa/tau_ampa +z/ms : siemens  
#dz/dt = -z/tau_ampa +g_syn*Tr_pre/ms : siemens
#g_syn=g_synpk/(tau_ampa/ms*exp(-1)) : siemens
#Tr_pre=.25*(tanh((t/ms-tspike/ms)/.005)-tanh((t/ms-(tspike/ms +transwdth/ms))/.005)) : 1
#spike : second
#''' 

synfireEqs = '''
dv/dt = ((((0.04)*v**2+(5)*v+140-u+I)+x2)*((1./taum)*ms))/ms : 1
du/dt = a5*(b5*v-u)/ms : 1
dx2/dt = ((-x2+y2/mV)*(1./taupsp)) : 1
dy2/dt = -y2*(1./taupsp)+25.27*mV/ms+(39.24*mV/ms**0.5)*xi : volt
I : 1
x : metre
y : metre
area : 1
previous_v : volt
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
v=msn_c
u +=msn_d
'''
reset3='''
v=c3
u +=d3
'''

reset4='''
v=c4
u+=d4
'''

synfireReset='''
v=c5
u+=d5
'''

#resetAmpa='''
#v=c5
#u+=d5
#'''

neuron_spacing1 = 50*umetre
width = n/4.0*neuron_spacing1

neuron_spacing2 = 500*umetre
width = n/4.0*neuron_spacing2

neuron_spacing3 = 20*umetre
width = n/4.0*neuron_spacing3

neuron_spacing4 = 20*umetre
width = n/4.0*neuron_spacing4

neuron_spacing5 = 20*umetre
width = n/4.0*neuron_spacing5

#%%
###############################################################################
############ neurongroups #####################################################
###############################################################################

# ActionSelection
ActionSelection=NeuronGroup(10,eqs4,threshold='v>30*mV',reset=reset4,method='euler')

ActionSelection.v = c4
ActionSelection.u = b4*c4
ActionSelection.x = 'i*neuron_spacing5'

# thalamo-cortical cells
n_groups = 10
group_size = 100
Cortex=NeuronGroup(N=n_groups*group_size, model=synfireEqs,threshold='v > 30',
                    reset=synfireReset, 
                    method='euler')
Cortex.v = c5
Cortex.u = b5*c5
Cortex.x = 'i*neuron_spacing1'
Cortex.y = 'i*neuron_spacing1'
Cortex.area = 100 # um, aribtrary

# motor cortex
MotorCortex=NeuronGroup(n, eqs,threshold='v > 30*mV',
                    reset=reset, 
                    method='euler')
MotorCortex.v = c
MotorCortex.u = b*c
MotorCortex.x = 'i*neuron_spacing1'
MotorCortex.y = 'i*neuron_spacing1'
MotorCortex.area = 100 # um, aribtrary

# DA neurons
DA = NeuronGroup(n,eqs,threshold='v>30*mV',reset=reset,method='euler')

DA.v = c
DA.u = b*c
DA.x = 'i*neuron_spacing1'

# medium spiny neurons 
D1=NeuronGroup(n,eqs2,threshold='v > 30*mV',
                    reset=reset2,
                    method='euler')
                    
D1.v = msn_c
D1.u = msn_b*msn_c
D1.x = 'i*neuron_spacing2'

D2=NeuronGroup(n,eqs2,threshold='v > 30*mV',
                    reset=reset2,
                    method='euler')
                    
D2.v = msn_c
D2.u = msn_b*msn_c
D2.x = 'i*neuron_spacing2'

#GPe neurons
GPe=NeuronGroup(n,eqs4,threshold='v>30*mV',reset=reset4,method='euler')

GPe.v = c4
GPe.u = b*c
GPe.x = 'i*neuron_spacing5'

# Hypothalamus
Hypothalamus=NeuronGroup(n,eqs4,threshold='v>30*mV',reset=reset4,method='euler')

Hypothalamus.v = c4
Hypothalamus.u = b4*c4
Hypothalamus.x = 'i*neuron_spacing5'

# SNr neurons
SNr=NeuronGroup(n,eqs3,threshold='v>30*mV', reset=reset3,method='euler')

SNr.v = c3
SNr.u = b3*c3
SNr.x = 'i*neuron_spacing3'

# STN neurons
STN=NeuronGroup(n,eqs4,threshold='v>20*mV',reset=reset4,method='euler') #ampaEqs,threshold='v>20*mV',reset=reset4,method='euler')

STN.v = c4#c5
STN.u = b4*c4#b5*c5/ms
STN.x = 'i*neuron_spacing5'

# thalamic neurons
Thalamus =NeuronGroup(n,eqs,threshold='v>30*mV',reset=reset,method='euler')

Thalamus.v = c
Thalamus.u = b*c
Thalamus.x = 'i*neuron_spacing4'


# poisson groups
#P = PoissonGroup(n, np.arange(100)*Hz + 10*Hz) # input to cortical neurons
P2 = PoissonGroup(n,np.arange(100)*Hz + 10*Hz) # input to STN neurons
P3 = PoissonGroup(n,np.arange(100)*Hz + 10*Hz) # input to thalamic neurons
P4 = PoissonGroup(n,np.arange(100)*Hz + 10*Hz) # input to SNr 
P5 = PoissonGroup(n,np.arange(100)*Hz + 10*Hz) # input to GPe

# Spike generator groups
#indices = array([0, 2, 1])
#times = array([1, 2, 3])*ms
#inp = SpikeGeneratorGroup(3, indices, times)

#indices = array(range(0,n))
#times = array(range(1,n+1))*ms
#inp = SpikeGeneratorGroup(n, indices, times)


#%%
####################################################################
############ synapses ##############################################
####################################################################

# STDP variables
taupre = taupost = 20*ms
wmax = 1 # weights can go from 0 to 1
Apre = 0.01
Apost = -Apre*taupre/taupost*1.05

# STDP eqs
stdp_eqs = '''
             w : 1
             dapre/dt = -apre/taupre : 1 (event-driven)
             dapost/dt = -apost/taupost : 1 (event-driven)
             '''
onpre='''
             v_post += w*mV
             apre += Apre
             w = clip(w+apost, 0, wmax)
             '''
onpreNeg='''
             v_post -= w*mV
             apre += Apre
             w = clip(w+apost, 0, wmax)
             '''      
onpreDA='''
             v_post *= w
             apre += Apre
             '''             
onpost='''
             apost += Apost
             w = clip(w+apre, 0, wmax)
             '''
onpostDA='''
             apost += Apost
         '''

#Asynapse = alpha_synapse(input='x',tau=10*ms,unit=amp)             

# Start weight maps based on distance
distWeights = '(exp(-(x_pre-x_post)**2/(2*width**2)))'

# Poisson input to cortex
#S2 = Synapses(P, Cortex, on_pre='v+=1')
#S2.connect(condition='i!=j',p=0.2)
#S2.delay = 2*ms

# Poisson input to STN; excitatory
S5 = Synapses(P2,STN,on_pre='v+=5*mV')
S5.connect(condition='i!=j',p=1)
S5.delay = 2*ms

# Poisson input to thalamus
S7 = Synapses(P3,Thalamus,on_pre='v+=1.75*mV')
S7.connect(condition='i!=j',p=1)
S7.delay = 2*ms

# Poisson input to SNr
P_SNr = Synapses(P4,SNr,on_pre='v+=1.75*mV')
P_SNr.connect(condition='i!=j',p=1)
P_SNr.delay = 2*ms

# Poisson input to GPe
P_GPe = Synapses(P5, GPe, on_pre='v+=1.75*mV')
P_GPe.connect(condition='i!=j',p=1)
P_GPe.delay = 2*ms

##############################################################################
# Spike generators
#feedforward_Cortex = Synapses(inp, Cortex,on_pre='v+=1.75')
#feedforward_Cortex.connect(j='i')

##############################################################################
# excitatory

# cortical to BG synapses; excitatory
Cortex_D1 = Synapses(Cortex,D1,stdp_eqs,on_pre=onpre,on_post=onpost)
Cortex_D1.connect(condition='i!=j',p=0.5)
Cortex_D1.delay = 10*ms
Cortex_D1.w = distWeights # input varies with distance

Cortex_D2 = Synapses(Cortex,D2,stdp_eqs,on_pre=onpre,on_post=onpost)
Cortex_D2.connect(condition='i!=j',p=0.5)
Cortex_D2.delay = 10*ms
Cortex_D2.w = distWeights # input varies with distance

# Excitatory projections from cortex to hypothalamus  
S_Cortex_Hypo = Synapses(Cortex,Hypothalamus,stdp_eqs,on_pre=onpre,on_post=onpost)      
S_Cortex_Hypo.connect(condition='i!=j',p=0.6)
S_Cortex_Hypo.delay = 5*ms
S_Cortex_Hypo.w= distWeights

# Excitatory projections from cortex to STN
S_Cortex_STN = Synapses(Cortex,STN,stdp_eqs,on_pre=onpre,on_post=onpost) #on_pre='''z_post += g_synpk''')#stdp_eqs,on_pre=onpre,on_post=onpost)
S_Cortex_STN.connect(condition='i!=j',p=0.5)
S_Cortex_STN.delay = 2.5*ms
#S_Cortex_STN.w= distWeights

# Excitatory Hypothalamic projections to DA
S_Hypo_DA = Synapses(Hypothalamus,DA,'w:1',on_pre='v=+10*mV')
S_Hypo_DA.connect(condition='i!=j',p=0.8)
S_Hypo_DA.delay = 5*ms
S_Hypo_DA.w=distWeights

# Excitatory STN to GPe
S_STN_GPe = Synapses(STN,GPe,stdp_eqs,on_pre=onpre,on_post=onpost)
S_STN_GPe.connect(condition='i!=j',p=0.2)
S_STN_GPe.delay = 2*ms
#S_STN_GPe.w = '(10*exp(-(x_pre-x_post)**2/(2*width**2)))' # input varies with distance

#STN to SNr; Excitatory
STN_SNr = Synapses(STN,SNr,stdp_eqs,on_pre=onpre,on_post=onpost)
STN_SNr.connect(condition='i!=j',p=0.2)
STN_SNr.delay = 1.5*ms
STN_SNr.w = distWeights # input varies with distance

# Thalamic input to cortex; excitatory 
Thal_Cortex = Synapses(Thalamus,Cortex,'w:1',on_pre='y2+=weight')
Thal_Cortex.connect(condition='i!=j',p=0.2)
Thal_Cortex.delay = 2*ms
Thal_Cortex.w = distWeights

#Cortex to ActionSelection; excitatory
MotorCortex_Action = Synapses(MotorCortex,ActionSelection,'w:1',on_pre='v+=2*mV')
MotorCortex_Action.connect(condition='i!=j',p=1)
MotorCortex_Action.delay = 5*ms
MotorCortex_Action.w = distWeights # input varies with distance

# Thalamus to Striatum; excitatory
Thal_D1 = Synapses(Thalamus,D1,stdp_eqs,on_pre=onpre,on_post=onpost)
Thal_D1.connect(condition='i!=j',p=0.3) # used lower probability than cortex
Thal_D1.delay = 10*ms # used same delay as cortex
Thal_D1.w = distWeights

Thal_D2 = Synapses(Thalamus,D2,stdp_eqs,on_pre=onpre,on_post=onpost)
Thal_D2.connect(condition='i!=j',p=0.3) # used lower probability than cortex
Thal_D2.delay = 10*ms # used same delay as cortex
Thal_D2.w = distWeights
###############################################################################3
# Inhibitory connections

# BG input to SNr; inhibitory
BG_SNr = Synapses(D1,SNr,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
BG_SNr.connect(condition='i!=j',p=1)
BG_SNr.delay = 4*ms
BG_SNr.w = distWeights # input varies with distance

# BG input to GPe; inhibitory
BG_GPe = Synapses(D2,GPe,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
BG_GPe.connect(condition='i!=j',p=0.34)
BG_GPe.delay = 5*ms
BG_GPe.w = distWeights # input varies with distance

# SNr to thalamic; Inhibitory
SNr_Thal = Synapses(SNr,Thalamus,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
SNr_Thal.connect(condition='i!=j',p=0.3)
SNr_Thal.delay = 4*ms
SNr_Thal.w = distWeights # input varies with distance

# GPe to STN; Inhibitory
GPe_STN = Synapses(GPe,STN,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
GPe_STN.connect(condition='i!=j',p=0.4)
GPe_STN.delay = 4*ms
GPe_STN.w = distWeights # input varies with distance

#Inhibitory projections to hypothalamus
#ActionSelection_Hypo = Synapses(ActionSelection,Hypothalamus,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
#ActionSelection_Hypo.connect(condition='i!=j',p=.7)
#ActionSelection_Hypo.delay = 5*ms
#ActionSelection_Hypo.w = distWeights # input varies with distance

##################################################################
# Modulatory 

# DA projections to D1
S_DA_D1 = Synapses(DA,D1,stdp_eqs,on_pre=onpreDA,on_post=onpostDA)
S_DA_D1.connect(condition='i!=j',p=0.1)
S_DA_D1.delay = 5*ms
S_DA_D1.w= distWeights

# DA projections to D2
S_DA_D2 = Synapses(DA,D2,stdp_eqs,on_pre=onpreDA,on_post=onpostDA)
S_DA_D2.connect(condition='i!=j',p=0.1)
S_DA_D2.delay = 5*ms
S_DA_D2.w= distWeights

################################################################
# Collaterals and recurrents

# Inhibitory SNr collaterals 
S_SNr_SNr = Synapses(SNr,SNr,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
S_SNr_SNr.connect(condition='i!=j',p=0.5)
S_SNr_SNr.delay = 1*ms
S_SNr_SNr.w = distWeights

# Inhibitory Striatal collaterals
D1_D1 = Synapses(D1,D1,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
D1_D1.connect(condition='i!=j',p=0.1)
D1_D1.delay = 1*ms

D2_D2 = Synapses(D2,D2,stdp_eqs,on_pre=onpreNeg,on_post=onpost)
D2_D2.connect(condition='i!=j',p=0.1)
D2_D2.delay = 1*ms

# Recurrent cortical connections 
#Cortex_Cortex = Synapses(Cortex,Cortex,stdp_eqs, on_pre=onpre,on_post=onpost)
#Cortex_Cortex.connect(condition='i!=j',p=0.1) # probability used in Laje & Buonomano, 2013
#Cortex_Cortex.delay = 1*ms

# Synfire cortical connections
Cortex_cortex = Synapses(Cortex,Cortex,on_pre='y2+=weight')
Cortex_cortex.connect(j='k for k in range((int(i/group_size)+1)*group_size, (int(i/group_size)+2)*group_size) '
'if i<N_pre-group_size')



#%%
####################################################################
############ LFP recorder ##########################################
####################################################################

#Cortex.run_regularly('previous_v = v', when='start')


#Ne = 5 # Number of electrodes
#sigma = 0.3 # Resistivity of extracellular field (0.3-0.4 S/m)
#lfp = NeuronGroup(Ne,model='''v : volt
#                              x : meter
#                              y : meter
#                              z : meter''')
#lfp.x = 7*mm 
#lfp.y = [1*mm, 2*mm, 4*mm, 8*mm, 16*mm]
## Synapses are normally executed after state update, so v-previous_v = dv
#LFPsynapses = Synapses(Cortex,lfp,model='''w : 1 #ohm*meter**2 (constant) # Weight in the LFP calculation
#                                 v_post = w*((1*metre-i*neuron_spacing1)*(v_pre-previous_v_pre)/dt) : volt (summed)''')
#LFPsynapses.connect(condition='i!=j',p=1)
#    # not quite area
#LFPsynapses.w = 'area_pre/(4*pi*sigma)/(((x_pre-x_post)/metre)**2+((y_pre-y_post)/metre)**2)' #(4*pi*sigma)
#LFPsynapses.summed_updaters['v_post'].when = 'after_groups'  # otherwise v and previous_v would be identical
#
#Mlfp = StateMonitor(lfp,'v',record=True)


#%%
####################################################################
############ monitors ##############################################
####################################################################

# neurons
ActionSpikes = SpikeMonitor(ActionSelection)
ActionTrace = StateMonitor(ActionSelection,('v'),record=True)

cortexSpikes = SpikeMonitor(Cortex)
cortexTrace = StateMonitor(Cortex, 'v', record=True)

D1_Spikes = SpikeMonitor(D1)
D1_Trace = StateMonitor(D1, ('v'), record=True)

D2_Spikes = SpikeMonitor(D2)
D2_Trace = StateMonitor(D2, ('v'), record=True)

DAspikes = SpikeMonitor(DA)
DAtrace = StateMonitor(DA,'v',record=True)

GPeSpikes = SpikeMonitor(GPe)
GPeTrace = StateMonitor(GPe, ('v'), record=True)

HypoSpikes  = SpikeMonitor(Hypothalamus)
HypoTrace = StateMonitor(Hypothalamus, ('v'), record=True)

SNrSpikes = SpikeMonitor(SNr)
SNrTrace = StateMonitor(SNr, ('v'), record=True)

STNSpikes = SpikeMonitor(STN)
STNTrace = StateMonitor(STN, ('v'), record=True)

thalamusSpikes = SpikeMonitor(Thalamus)
thalamusTrace = StateMonitor(Thalamus, ('v'), record=True)

#%%
####################################################################
############# Network Operations ###################################
####################################################################

# This will cause LTP at neurons that have had DA and cortical input
Cortex_D1mon = StateMonitor(Cortex_D1,'w', record=Cortex_D1['i!=j'])
Cortex_D2mon = StateMonitor(Cortex_D2,'w', record=Cortex_D2['i!=j'])

#%%
#####333
# temp seperation... think about making this not the inverse, but the same STDP as LTP, except divide by DA signal instead of multiply

@network_operation(dt=10*ms)
def DA_LTD():
    CortexIndex = cortexSpikes.i[cortexSpikes.t > defaultclock.t-window]  # index for the active neuron
    DAindex = DAspikes.i[DAspikes.t > defaultclock.t-window]
    bothIndex = intersect1d(array(CortexIndex),array(DAindex))
    synapseIndexArray = Cortex_D2.i[0:np.size(Cortex_D2.i)]
    synapseIndexMat = []
    for j in range(0,np.size(bothIndex)):
        if np.size(bothIndex) > 0:
            synapseIndex = np.where(synapseIndexArray==bothIndex[j])[0]
            synapseIndexMat.extend(synapseIndex)    
    #print synapseIndexMat
    #Cortex_D2.w_[synapseIndexMat] *= DAltdRate
    Cortex_D2.w_[synapseIndexMat] = clip(Cortex_D2.w_[synapseIndexMat]+(-Cortex_D2.apost_[synapseIndexMat]*DAltdRate),0,wmax)
    
  
@network_operation(dt=10*ms)
def DA_LTP():
    CortexIndex = cortexSpikes.i[cortexSpikes.t > defaultclock.t-window]  # index for the active neuron
    DAindex = DAspikes.i[DAspikes.t > defaultclock.t-window]
    bothIndex = intersect1d(array(CortexIndex),array(DAindex))
    synapseIndexArray = Cortex_D1.i[0:np.size(Cortex_D1.i)]
    synapseIndexMat = []
    for j in range(0,np.size(bothIndex)):
        if np.size(bothIndex) > 0:
            synapseIndex = np.where(synapseIndexArray==bothIndex[j])[0]
            synapseIndexMat.extend(synapseIndex)    
    #print synapseIndexMat
    #Cortex_D1.w_[synapseIndexMat] *= DAltpRate #1/(1+math.exp(DAltpRate))                  
    Cortex_D1.w_[synapseIndexMat] = clip(Cortex_D1.w_[synapseIndexMat]+(Cortex_D1.apost_[synapseIndexMat]*DAltpRate),0,wmax)       
    #Cortex_D1.w_[synapseIndexMat] = clip(Cortex_D1.w_[synapseIndexMat]*DAltpRate),0,wmax)       
    
    

# write an operation that will decrease current as homeostasis is reached
@network_operation(dt=10*ms)
def motivation():
    Hypothalamus.I *= (1000-ActionSpikes.count[0])/1000 #clip(-ActionSpikes.count[0]*feedbackConstant,0,1000*nA) # decreases based upon internal state
    # consider changing this to an exponential decay function

#%%
####################################################################
############ run ###################################################
####################################################################

Hypothalamus.I = 100*nA # motivation input
Thalamus.I = 50*nA#'12*nA + 12*i/n*nA' # some salient sensory input
#inp.set_spikes(cortexSpikes.i, cortexSpikes.t + runtime) # feed previous output of Cortex as input to group
#DA.I = 100*nA
#Cortex.I = 100*nA
run(duration, report='text')
Hypothalamus.I = 100*nA # motivation input
Thalamus.I = 0*nA#'12*nA + 12*i/n*nA' # some salient sensory input
#inp.set_spikes(cortexSpikes.i, cortexSpikes.t + runtime) # feed previous output of Cortex as input to group
run(duration, report='text')
Hypothalamus.I = 100*nA # motivation input
Thalamus.I = 50*nA#'12*nA + 12*i/n*nA' # some salient sensory input
#inp.set_spikes(cortexSpikes.i, cortexSpikes.t + runtime) # feed previous output of Cortex as input to group
run(duration, report='text')
#Hypothalamus.I = 100*nA # motivation input
#Thalamus.I = 0*nA#'12*nA + 12*i/n*nA' # some salient sensory input
#inp.set_spikes(cortexSpikes.i, cortexSpikes.t + runtime) # feed previous output of Cortex as input to group
#run(duration, report='text')
#Hypothalamus.I = 100*nA # motivation input
#Thalamus.I = 50*nA#'12*nA + 12*i/n*nA' # some salient sensory input
#inp.set_spikes(cortexSpikes.i, cortexSpikes.t + runtime) # feed previous output of Cortex as input to group
#run(duration, report='text')
#Hypothalamus.I = 100*nA # motivation input
#Thalamus.I = 0*nA#'12*nA + 12*i/n*nA' # some salient sensory input
#inp.set_spikes(cortexSpikes.i, cortexSpikes.t + runtime) # feed previous output of Cortex as input to group
#run(duration, report='text')


#%%
####################################################################
############ plotting ##############################################
####################################################################
figure(1)
plot(cortexTrace.t/ms,cortexTrace.v[1]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('Cortical Voltage')

figure(2)
plot(D1_Trace.t/ms,D1_Trace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('D1_MSN Voltage')

figure(3)
plot(SNrTrace.t/ms,SNrTrace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('SNr Voltage')

figure(4)
plot(thalamusTrace.t/ms,thalamusTrace.v[70]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('Thalamic Voltage')

figure(5)
plot(GPeTrace.t/ms,GPeTrace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('GPe Voltage')

figure(6)
plot(STNTrace.t/ms,STNTrace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('STN Voltage')

figure(7)
plot(DAtrace.t/ms,DAtrace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('DA Voltage')

figure(8)
plot(HypoTrace.t/ms,HypoTrace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('Hypothalamus Voltage')

figure(9)
plot(ActionTrace.t/ms,ActionTrace.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('Brainstem Voltage')

figure(10)
plot(Cortex_D1mon.t/ms, Cortex_D1mon.w[0])
xlabel('Time(ms)')
ylabel('Weights')
title('Cortical D1 type weight')

figure(11)
plot(Cortex_D1mon.t/ms, Cortex_D2mon.w[0])
xlabel('Time(ms)')
ylabel('Weights')
title('Cortical D2 type weight')

figure(12)
plot(cortexSpikes.t/ms, 1.0*cortexSpikes.i/group_size, '.')
plot([0, duration/ms], np.arange(n_groups).repeat(2).reshape(-1, 2).T, 'k-')
ylabel('group number')
yticks(np.arange(n_groups))
xlabel('time (ms)')
show()


#figure(12)
#subplot(211)
#for i in range(10):
#    plot(cortexTrace.t/ms,cortexTrace[i].v/mV)
#ylabel('$V_m$ (mV)')
#subplot(212)
#for i in range(5):
#    plot(cortexTrace.t/ms,Mlfp.v[i]/mV)
#ylabel('LFP (mV)')
#xlabel('Time (ms)')
#show()

