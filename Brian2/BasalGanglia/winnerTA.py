# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 10:05:50 2016

@author: elijah
"""

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

eqs2 = '''
dv/dt = ((0.04)*v**2+(5)*v+140-u+I)/ms : 1
du/dt = a*(b*v-u)/ms : 1
I : 1
'''

reset='''
v=c
u += d
'''

#%%
# Nuerons
# Actions
LeverPress = NeuronGroup(10,eqs2,threshold='v>20',reset=reset,method='euler')
WrongLeverPress = NeuronGroup(10,eqs2,threshold='v>20',reset=reset,method='euler')
NoLeverPress = NeuronGroup(10,eqs2,threshold='v>20',reset=reset,method='euler')


#%% Synapses
# Action to Action
# this is building a winner take all circuit
LeverPress_NoPress = Synapses(LeverPress,NoLeverPress,on_pre='v-=10') # These are inhibitory with longer delays
LeverPress_NoPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_NoPress.delay = 5*ms

LeverPress_NoPressEx = Synapses(LeverPress,NoLeverPress,on_pre='v+=5') # excitatory with lower delays 
LeverPress_NoPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_NoPressEx.delay = 2*ms

LeverPress_WrongPress = Synapses(LeverPress,WrongLeverPress,on_pre='v-=10')
LeverPress_WrongPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_WrongPress.delay = 5*ms

LeverPress_WrongPressEx = Synapses(LeverPress,WrongLeverPress,on_pre='v+=5')
LeverPress_WrongPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_WrongPressEx.delay = 2*ms

# this is building a winner take all circuit
WrongPress_NoPress = Synapses(WrongLeverPress,NoLeverPress,on_pre='v-=10')
WrongPress_NoPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_NoPress.delay = 5*ms

WrongPress_NoPressEx = Synapses(WrongLeverPress,NoLeverPress,on_pre='v+=5')
WrongPress_NoPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_NoPressEx.delay = 2*ms

WrongPress_Press = Synapses(WrongLeverPress,LeverPress,on_pre='v-=10')
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

#%% Record
ActionSpikes = SpikeMonitor(LeverPress)
NoActionSpikes = SpikeMonitor(NoLeverPress)
WrongActionSpikes = SpikeMonitor(WrongLeverPress)

#%% Run
NoLeverPress.I = 0
LeverPress.I = 100
WrongLeverPress.I = 0
run(duration,report='text')

#%% Plot
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

figure()
plot(ActionSpikes.t/ms,ActionSpikes.i,'xb')
plot(NoActionSpikes.t/ms,NoActionSpikes.i,'or')
plot(WrongActionSpikes.t/ms,WrongActionSpikes.i,'.g')
xlabel('Time(ms)')
ylabel('Neuron Number')
title('Brainstem Spikes')
    