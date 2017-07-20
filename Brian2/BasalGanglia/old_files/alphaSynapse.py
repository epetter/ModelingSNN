# -*- coding: utf-8 -*-
"""
Created on Thu Oct 06 08:49:48 2016

@author: Elijah
"""
from brian2 import *
start_scope()
duration = 100*ms
vt = 30*mV         # Spiking threshold
memc = 200.0*pfarad  # Membrane capacitance
bgcurrent = 200*pA   # External current
tau_m=6*ms
tau_ampa=8*ms
g_synpk=0.01*nsiemens
transwdth=1.0*ms

#%%
#################################################################
############### Equations #######################################
#################################################################

eqs_neurons='''
dv/dt= -v/tau_m +  I*volt/second : volt (unless refractory)
I : 1
'''

eqs_neurons2='''
dv/dt= -v/tau_m + g_ampa*volt/memc  : volt (unless refractory)
dg_ampa/dt = -g_ampa/tau_ampa +z/ms : siemens
dz/dt = -z/tau_ampa +g_syn*Tr/ms : siemens
g_syn=g_synpk/(tau_ampa/ms*exp(-1)) : siemens
Tr=.25*(tanh((t/ms-tspike/ms)/.005)-tanh((t/ms-(tspike/ms +transwdth/ms))/.005)) : 1

I : 1
tspike:second
'''
#%%
#################################################################
############### Neurons #########################################
#################################################################


N = NeuronGroup(1, model=eqs_neurons, threshold='v > vt',
                      reset='v=-60*mV',refractory='25*ms',
                      method='euler')

N.I=17
N.v=0*mV
#neurons.tspike=-100.0*ms  #needed so a spike does not happen at time 0
# Comment these two lines out to see what happens without Synapses

N2=NeuronGroup(1,model=eqs_neurons2, threshold='v > vt',
                      reset='v=-60*mV;g_ampa=0*nS;tspike=t',refractory='2*ms',
                      method='euler')

N2.v=0.0*mV
N2.tspike=-100.0*ms

P1 = PoissonGroup(100,np.arange(100)*Hz + 10*Hz)


#%%
#################################################################
############### Synapses ########################################
#################################################################
S=Synapses(N,N2, on_pre='''
 z_post += g_synpk
 ''')
S.connect(i=0,j=0)

#S2=Synapses(P1,N,on_pre='v+=125*mV')
#S2.connect(i=0,j=0)

#%%
#################################################################
############### Monitors #########################################
#################################################################
M = StateMonitor(N, ('v'), record=True)
sm = SpikeMonitor(N)
M2 = StateMonitor(N2, ('v','g_ampa','Tr'), record=True)
sm2 = SpikeMonitor(N2)

#%%
#################################################################
############### Run #############################################
#################################################################
run(duration)

#%%
#################################################################
############### Plotz ###########################################
#################################################################
figure(1)
title('alpha synapse')
subplot(3,1,1)
plot(M.t/ms, M.v[0]/mV, '-b')
xlim(0,duration/ms)
subplot(3,1,2)
plot(M2.t/ms,M2.g_ampa[0]/nsiemens,'-g', lw=2)
subplot(3,1,3)
plot(M2.t/ms,M2.v[0]/mV,'-b')
show()