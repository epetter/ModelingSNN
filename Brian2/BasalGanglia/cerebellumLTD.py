# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 15:47:32 2016

@author: Elijah
"""

from pylab import *
from brian2 import *

start_scope()

#neuron_spacing1 = 50*umetre
#width = N_PKJ/4.0*neuron_spacing1

from pylab import *
from brian2 import *

class LTD(SpikeMonitor):
    '''
    Implements LTD on PF-PKJ synapses driven by CF activity
    '''
    def __init__(self, IO, GR, PF_PKJ, CF_PKJ, ltd_decay=.995, window=50*ms):
        '''
        IO: NeuronGroup of inferior olive neurons
        GR: NeuronGroup of granule cells
        PF_PKJ: Synapses object of PF-PKJ synapses
        CF_PKJ: Synapses object of CF-PKJ synapses
        ltd_decay: the multiplicative decay constant to decay PF-PKJ
                   synapses by.
        window: the time window to depress PF-PKJ synapses for GR spikes 
                that precede IO spikes.
        '''
        self.IO = IO
        self.GR = GR
        self.PF_PKJ = PF_PKJ
        self.CF_PKJ = CF_PKJ
        self.ltd_decay = ltd_decay
        self.ltd_bins = int(window/defaultclock.dt)
        SpikeMonitor.__init__(self, IO)
        
    def propagate(self, spikes):
        '''
        Depress PF-PKJ synapses for GRs that spiked that synapse on
        PKJs that received CF spike.
        '''
        if len(spikes):
            # PKJs that received CF spike
            pkj_inds = LTD.postsynaptic_indexes(spikes,self.CF_PKJ)
            
            # GRs that synapse on pkj_inds
            gr_inds = LTD.presynaptic_indexes(pkj_inds,self.PF_PKJ)
            
            # GRs that spiked in the past ltd_bins
            gr_spiked_inds = self.GR.LS[:self.ltd_bins]

            # GRs that spiked and synapse on PKJs
            gr_ltd_inds = intersect1d(array(gr_inds),array(gr_spiked_inds))
            if not len(gr_ltd_inds):
                return
            
            # modify PF_PKJ synaptic strength
            self.PF_PKJ.w[gr_ltd_inds,pkj_inds] *= self.ltd_decay

    @staticmethod
    def presynaptic_indexes(neuron_inds,synapses):
        '''
        Map from target neuron indexes to source neuron indexes
        in synapses
        '''
        synapse_indexes = hstack(array(synapses.i)[neuron_inds])
        return synapses.presynaptic[synapse_indexes]
    
    @staticmethod
    def postsynaptic_indexes(neuron_inds,synapses):
        '''
        Map from source neuron indexes to target neuron indexes
        in synapses
        '''
        synapse_indexes = hstack(array(synapses.i)[neuron_inds])
        return synapses.postsynaptic[synapse_indexes]


PF_PKJ.i
# Number of neurons
N_GR = 10000
N_PKJ = 16
N_IO  = 80 

# variables
# thalamo-cortical Iz
c = -66*mV
b= 0.2/ms
a = 0.02/ms
d = 8*mV/ms

# eqs
eqs =  '''
dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u+I/nF : volt
du/dt = a*(b*v-u) : volt/second
I : amp
x : metre
'''

reset='''
v=c
u += d
'''

# Neuron groups
GR = NeuronGroup(N_GR,eqs,threshold='v>30*mV',reset=reset,method='euler') # Granule cells -- their axons are PFs
PKJ = NeuronGroup(N_PKJ,eqs,threshold='v>30*mV',reset=reset,method='euler') #
#PKJ.x = 'i*neuron_spacing1'
IO = NeuronGroup(N_IO,eqs,threshold='v>30*mV',reset=reset,method='euler') # Inferior Olive neurons -- their axons are CFs 

# Synapses
CF_PKJ = Synapses(IO, PKJ,'w:1',on_pre='v+=w*mV')
CF_PKJ.connect(condition='i!=j',p=1)
CF_PKJ.w = 0.9


PF_PKJ = Synapses(GR,PKJ,'w:1',on_pre='v+=w*mV')
PF_PKJ.connect(condition='i!=j',p=1)# GR -> PKJ randomly, sparsely
PF_PKJ.w = 0.3 #'(exp(-(x_pre-x_post)**2/(2*width**2)))' # input varies with distance#monitor = LTD(IO,GR,PF_PKJ,CF_PKJ)
monitor = LTD(IO,GR,PF_PKJ,CF_PKJ)
spikes = SpikeMonitor(PKJ)
pkj_inds = LTD.postsynaptic_indexes(spikes.i,CF_PKJ)
PF_PKJ.w = LTD.propagate(monitor,spikes)

#PF_PKJ = monitor.PF_PKJ

#monitor.propagate(monitor)
# connect IO -> PKJ, each CF synapses on 2 PKJ; no PKJ receives more than one CF

S = StateMonitor(IO, 'v', record=True)
S2 = StateMonitor(PKJ, 'v', record=True)
S3 = StateMonitor(GR, 'v', record=True)
spikez = SpikeMonitor(PKJ)

IO.I = 100*nA
GR.I = 100*nA
run(400*ms, report='text')

figure(1)
plot(S2.t/ms,S2.v[0]/mV)
xlabel('Time(ms)')
ylabel('Voltage(mV)')
title('PKJ Voltage')

figure(2)
plot(S.t/ms, S.v[0]/mV)
title('IO voltage')

figure(3)
plot(S3.t/ms, S3.v[0]/mV)
title('GR voltage')

figure(4)
plot(PF_PKJ.w)

figure(5)
plot(monitor.PF_PKJ.w)