# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 18:58:29 2017

@author: Elijah
"""

import numpy as np

def spike_covariance(**kwargs):
    # this function will calculat covariance between n firing rates
        # example
    #spike_cov, firing_rates = spike_covariance(f1=SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,
    #                             f2=D1pop.smooth_rate(window='gaussian',width=binSize)/Hz)
    return np.cov(np.array([x for key, x in kwargs.iteritems()])), np.array([x for key, x in kwargs.iteritems()])
    
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
    
def visualize_stdp(tau_pre,A_pre,A_post):     
    tau_post = tau_pre 
    delta_t = linspace(-500, 500, 1000)*ms
    W = where(delta_t<0, A_pre*exp(delta_t/tau_pre), A_post*exp(-delta_t/tau_post))
    plot(delta_t/ms, W)
    xlabel(r'$\Delta t$ (ms)')
    ylabel('W')
    ylim(-A_post, A_post)
    axhline(0, ls='-', c='k')    
        