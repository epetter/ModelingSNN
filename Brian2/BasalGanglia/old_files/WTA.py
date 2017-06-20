# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 16:29:37 2016

@author: Elijah
"""
from brian2 import *
import numpy as np
start_scope()
duration = 1000*ms


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

reset='''
v=c
u += d
'''

reset3='''
v=c2
u += d2
'''

LeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
NoLeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
WrongLeverPress = NeuronGroup(2,eqs2,threshold='v>30',reset=reset,method='euler')
#InhInter = NeuronGroup(5,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterNL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')
InhInterWL = NeuronGroup(2,eqs3,threshold='v>20',reset=reset3,method='euler')

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

LeverPress.I = 50
ActionSpikes = SpikeMonitor(LeverPress)
NoActionSpikes = SpikeMonitor(NoLeverPress)
WrongActionSpikes = SpikeMonitor(WrongLeverPress)
InhSpikes = SpikeMonitor(InhInterL)
run(duration,report='text')
LeverPress.I = 10
NoLeverPress.I = 50
run(duration,report='text')
LeverPress.I = 10
WrongLeverPress.I = 50
run(duration,report='text')

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

Lfr,Lbin = calculate_FR(ActionSpikes,binSize=100*ms,timeWin=duration*3/second)
NLfr,NLbin = calculate_FR(NoActionSpikes,binSize=100*ms,timeWin=duration*3/second)
WLfr,WLbin = calculate_FR(WrongActionSpikes,binSize=100*ms,timeWin=duration*3/second)
figure()
plot(range(0,int(duration*3/ms)-360,100),Lbin,'r')
plot(range(0,int(duration*3/ms)-360,100),NLbin,'b')
plot(range(0,int(duration*3/ms)-360,100),WLbin,'g')
xlabel('Time(ms)')