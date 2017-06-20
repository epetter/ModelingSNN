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
LeverPress = NeuronGroup(2,eqs2,threshold='v>20',reset=reset,method='euler')
NoLeverPress = NeuronGroup(2,eqs2,threshold='v>20',reset=reset,method='euler')
WrongLeverPress = NeuronGroup(2,eqs2,threshold='v>20',reset=reset,method='euler')
InhInter = NeuronGroup(1,eqs2,threshold='v>20',reset=reset,method='euler')

Cortex = NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')

Cortex.v = c
Cortex.u = b*c


D1 = NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')

D1.v = c
D1.u = b*c

SNr = NeuronGroup(n,eqs2,threshold='v>30',reset=reset,method='euler')

SNr.v = c
SNr.u = b*c


#%% Synapses
# Action to Action
# this is building a winner take all circuit
LeverPress_NoPressEx = Synapses(LeverPress,NoLeverPress,on_pre='v+=5') # excitatory with lower delays 
LeverPress_NoPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_NoPressEx.delay = 5*ms

LeverPress_WrongPressEx = Synapses(LeverPress,WrongLeverPress,on_pre='v+=5')
LeverPress_WrongPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_WrongPressEx.delay = 5*ms

LeverPress_PressEx = Synapses(LeverPress,LeverPress,on_pre='v+=5') # excitatory with lower delays 
LeverPress_PressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_PressEx.delay = 5*ms

LeverPress_Inh = Synapses(LeverPress,InhInter,on_pre='v+=5')
LeverPress_Inh.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
LeverPress_Inh.delay = 5*ms


# Wrong lever Press
WrongPress_NoPressEx = Synapses(WrongLeverPress,NoLeverPress,on_pre='v+=5')
WrongPress_NoPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_NoPressEx.delay = 5*ms

WrongPress_PressEx = Synapses(WrongLeverPress,LeverPress,on_pre='v+=5')
WrongPress_PressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_PressEx.delay = 5*ms

WrongPress_WrongPressEx = Synapses(WrongLeverPress,WrongLeverPress,on_pre='v+=5')
WrongPress_WrongPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_WrongPressEx.delay = 5*ms

WrongPress_Inh = Synapses(WrongLeverPress,InhInter,on_pre='v+=5')
WrongPress_Inh.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#LeverPress_NoPress.connect(condition='i!=j',p=1,skip_if_invalid=True)
WrongPress_Inh.delay = 5*ms

# No Lever Press 
NoPress_PressEx = Synapses(NoLeverPress,LeverPress,on_pre='v+=5')
NoPress_PressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#NoPress_Press.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_PressEx.delay = 5*ms

NoPress_WrongPressEx = Synapses(NoLeverPress,WrongLeverPress,on_pre='v+=5')
NoPress_WrongPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#NoPress_Press.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_WrongPressEx.delay = 5*ms

NoPress_NoPressEx = Synapses(NoLeverPress,NoLeverPress,on_pre='v+=5')
NoPress_NoPressEx.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
#NoPress_Press.connect(condition='i!=j',p=1,skip_if_invalid=True)
NoPress_NoPressEx.delay = 5*ms

NoPress_Inh = Synapses(NoLeverPress,InhInter,on_pre='v+=5')
NoPress_Inh.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
NoPress_Inh.delay = 5*ms


# Inhibitory Interneuron
InhPress = Synapses(InhInter,LeverPress,on_pre='v=-6')
InhPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
InhPress.delay = 2*ms

InhWrongPress = Synapses(InhInter,WrongLeverPress,on_pre='v=-6')
InhWrongPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
InhWrongPress.delay = 2*ms

InhNoPress = Synapses(InhInter,NoLeverPress,on_pre='v=-6')
InhNoPress.connect(j='k for k in range(i-w2, i+w2) if rand()<1', skip_if_invalid=True) 
InhNoPress.delay = 2*ms


#%% Cortex To WTA circuit
Cortex_Lever = Synapses(Cortex,LeverPress,on_pre='v+=10')
Cortex_Lever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
Cortex_Lever.delay = 25*ms
Cortex_NoLever = Synapses(Cortex,NoLeverPress,on_pre='v+=10')
Cortex_NoLever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
Cortex_NoLever.delay = 25*ms
Cortex_WrongLever = Synapses(Cortex,WrongLeverPress,on_pre='v+=10')
Cortex_WrongLever.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True) 
Cortex_WrongLever.delay = 25*ms

# Cortex D1
Cortex_D1 = Synapses(Cortex,D1,on_pre='v=+10')
Cortex_D1.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True)
Cortex_D1.delay = 10*ms

D1_SNr = Synapses(D1,SNr,on_pre='v-=10')
D1_SNr.connect(j='k for k in range(i-w2, i+w2) if rand()<0.5', skip_if_invalid=True)
D1_SNr.delay = 4*ms

SNr_Action =  Synapses(SNr,LeverPress,on_pre='v-=10')
SNr_Action.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True)
SNr_Action.delay = 5*ms

#SNr_NoAction = Synapses(SNr,NoLeverPress,on_pre='v-=10')
#SNr_NoAction.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True)
#SNr_NoAction.delay = 5*ms

#SNr_WrongAction = Synapses(SNr,WrongLeverPress,on_pre='v-=10')
#SNr_WrongAction.connect(j='k for k in range(i-w2, i+w2) if rand()<0.2', skip_if_invalid=True)
#SNr_WrongAction.delay = 5*ms

#%% Record
ActionSpikes = SpikeMonitor(LeverPress)
NoActionSpikes = SpikeMonitor(NoLeverPress)
WrongActionSpikes = SpikeMonitor(WrongLeverPress)

#%% Run
#NoLeverPress.I = 0
#LeverPress.I = 10
#WrongLeverPress.I = 0
#run(duration,report='text')
#NoLeverPress.I = 0
#LeverPress.I = 0
#WrongLeverPress.I = 0
#run(duration,report='text')
#NoLeverPress.I = 5
#LeverPress.I = 5
#WrongLeverPress.I = 5
#run(duration,report='text')
#NoLeverPress.I = 0
#LeverPress.I = 0
#WrongLeverPress.I = 0
#run(duration,report='text')
#NoLeverPress.I = 5
#LeverPress.I = 5
#WrongLeverPress.I = 10
Cortex.I = 150
run(duration,report='text')


print len(NoActionSpikes.t)
print len(ActionSpikes.t)
print len(WrongActionSpikes.t)
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
ylim(-1,2)
xlabel('Time(ms)')
ylabel('Neuron Number')
title('Brainstem Spikes')
    