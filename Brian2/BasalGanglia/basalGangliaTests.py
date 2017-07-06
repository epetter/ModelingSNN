# -*- coding: utf-8 -*-
"""
Created on Wed Jul 05 11:11:34 2017

@author: elijah
"""

from brian2 import *
import numpy as np
import scipy as sp
import scipy.stats
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'C:\Users\elijah\Documents\GitHub\ModelingSNN\Brian2\BasalGanglia')
from basalGangliaNetwork import *

#%% Run and analyze

def sequence():  # reproduce figure 3 in humphries et al., 2006        
   # population activity

   CortexPopL = PopulationRateMonitor(CortexL)
   CortexPopNL = PopulationRateMonitor(CortexNL)
   CortexPopWL = PopulationRateMonitor(CortexWL)
   D1popL = PopulationRateMonitor(D1_L)
   D1popNL = PopulationRateMonitor(D1_NL)
   D1popWL = PopulationRateMonitor(D1_WL)
   SNrPopL = PopulationRateMonitor(SNrL)
   SNrPopNL = PopulationRateMonitor(SNrNL)
   SNrPopWL = PopulationRateMonitor(SNrWL)
   
   ta_CortexL = TimedArray([5,10,10],dt=sequence_duration)
   ta_CortexNL = TimedArray([5,5,40],dt=sequence_duration)
   ta_CortexWL = TimedArray([5,5,5],dt=sequence_duration)
  
   @network_operation(dt=sequence_duration)
   def update_sequence(t):
       CortexL.I = ta_CortexL(t)
       CortexNL.I = ta_CortexNL(t)
       CortexWL.I = ta_CortexWL(t)
   run(sequence_duration*3,report='text',report_period=report_time)
   
   figure()
   plot(SNrPopL.t/ms,SNrPopL.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(SNrPopNL.t/ms,SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(SNrPopWL.t/ms,SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('SNr Firing Rates')
   legend('R2U')
   plt.savefig(save_root + 'SequenceTest_SNrfiringRate.png')
    
   figure()
   plot(CortexPopL.t/ms,CortexPopL.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(CortexPopNL.t/ms,CortexPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(CortexPopWL.t/ms,CortexPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('Cortex Firing Rates')
   legend('R2U')
   plt.savefig(save_root + 'SequenceTest_CortexfiringRate.png')
   
   figure() 
   plot(D1popL.t/ms,D1popL.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(D1popNL.t/ms,D1popNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(D1popWL.t/ms,D1popWL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('D1 pop firing rates')
   legend('R2U')
   plt.savefig(save_root + 'SequenceTest_D1firingRate.png')

def pop_firing(): # reproduce figure 2 in Humphries et al., 2006   
   # population activity
   GPePop = PopulationRateMonitor(GPe_L)
   SNrPop = PopulationRateMonitor(SNrL)
   STNpop = PopulationRateMonitor(STN)
   
   run (pop_duration,report='text',report_period=report_time) 

   STN_FR = mean(STNpop.smooth_rate(window='gaussian',width=binSize)/Hz)
   GPe_FR = mean(GPePop.smooth_rate(window='gaussian',width=binSize)/Hz)
   SNr_FR = mean(SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz)
   
   STNemp = 10;
   GPeEmp = 29;
   SNrEmp = 26

   print 'STN'   
   print STN_FR
   print 'GPe'
   print GPe_FR
   print 'SNr'
   print SNr_FR

   figure()
   #plt.legend('STN','GPe','SNr')
   plt.plot(1,STN_FR,'ob')
   plt.plot(2,GPe_FR,'ob')
   plt.plot(3,SNr_FR,'ob')
   plt.plot(1,STNemp,'g+',mew=10, ms=20)
   plt.plot(2,GPeEmp,'g+', mew=10, ms=20)
   plt.plot(3,SNrEmp,'g+', mew=10, ms=20)
   plt.errorbar(1,STNemp,(14-8)/2,capsize=11, elinewidth=2)
   plt.errorbar(2,GPeEmp,(36-26)/2,capsize=11, elinewidth=2)
   plt.errorbar(3,SNrEmp,(32-21)/2,capsize=11, elinewidth=2)
   xlim(0,4)
   ylim(8,40)
   title('Mean Tonic Firing Rates')
   plt.savefig(save_root + 'PopTest_NucleifiringRate.png')
   
   print 'Population Firing Results'
   
   if (STN_FR > 8) & (STN_FR < 14):
      print 'STN firing rate passed test'
   else:
      print 'ERROR!! STN firing rate not reasonable'
      
   if (GPe_FR > 26) & (GPe_FR < 36):
       print 'GPe firing rate passed test'
   else:
      print 'ERROR!! GPe firing rate not reasonable'       

   if (SNr_FR > 21) & (SNr_FR < 32):
       print 'SNr firing rate passed test'
   else:
      print 'ERROR!! SNr firing rate not reasonable'           
      
def cortex_D1_action():
    CortexL_D1L.w = CortexD1_start * 2 # overwriting starting weights 
    
    # population monitors 
    CortexPop = PopulationRateMonitor(CortexL)
    D1Lpop = PopulationRateMonitor(D1_L)
    D1NLpop = PopulationRateMonitor(D1_NL)
    D1WLpop = PopulationRateMonitor(D1_WL) 
    SNrPop = PopulationRateMonitor(SNrL)
    SNrPopNL = PopulationRateMonitor(SNrNL)
    SNrPopWL = PopulationRateMonitor(SNrWL)
    
    CortexL.I = 5
    CortexNL.I = 5
    CortexWL.I = 5
    run(cortex_D1_duration,report='text',report_period=report_time)
    
    # smoothed rates 
    SNrL_smoothFR = SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz
    SNrNL_smoothFR =SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz
    SNrWL_smoothFR =SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz
        
    figure()
    plot(SNrPop.t/ms,SNrL_smoothFR,'r')
    plot(SNrPopNL.t/ms,SNrNL_smoothFR,'b')
    plot(SNrPopWL.t/ms,SNrWL_smoothFR,'g')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('SNr Firing Rates')
    legend('RUU')
    plt.savefig(save_root + 'CortexD1_SNrfiringRate.png')
    
    figure()
    plot(CortexPop.t/ms,CortexPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('Cortical Firing Rates')
    legend('R')
    plt.savefig(save_root + 'CortexD1_CortexfiringRate.png')
    
    figure()
    plot(D1Lpop.t/ms,D1Lpop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
    plot(D1NLpop.t/ms,D1NLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
    plot(D1WLpop.t/ms,D1WLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
    xlabel('Time(ms)')
    ylabel('Firing Rate')
    title('D1 Firing Rates')
    legend('R2U')   
    plt.savefig(save_root + 'CortexD1_D1firingRate.png')
    
    avg_D1L = np.mean(CortexL_D1L.w)
    avg_D1NL = np.mean(CortexNL_D1NL.w)
    avg_D1WL =  np.mean(CortexWL_D1WL.w)

    L_actions = np.sum(SNrPop.rate < SNr_thresh)
    NL_actions = np.sum(SNrPopNL.rate < SNr_thresh)
    WL_actions = np.sum(SNrPopWL.rate < SNr_thresh)
    
    actions = [L_actions, NL_actions, WL_actions]
    
    avg_FR = [np.mean(SNrL_smoothFR),np.mean(SNrNL_smoothFR),np.mean(SNrWL_smoothFR)]
    if np.where(avg_FR == np.min(avg_FR))[0] == 0:
       print 'cortex_D1_action test passed!!'
    else:
       print 'cortex_D1_action test FAILED!!'
       
    return actions 
    
def learn_action(): 
   # independent E/I Poisson inputs
   p1 = PoissonInput(CortexL, 'v', N=20, rate=2*Hz, weight=100)
   p2 = PoissonInput(CortexNL, 'v', N=20, rate=2*Hz, weight=100)
   p3 = PoissonInput(CortexWL, 'v', N=20, rate=2*Hz, weight=100)
   
   #%% Record
   SNrLspikes = SpikeMonitor(SNrL)
   SNrNLspikes = SpikeMonitor(SNrNL)
   SNrWLspikes = SpikeMonitor(SNrWL)
   STNspikes = SpikeMonitor(STN)
   CortexLspikes= SpikeMonitor(CortexL)
   CortexNLspikes = SpikeMonitor(CortexNL)
   CortexWLspikes = SpikeMonitor(CortexWL)
   DASpikes = SpikeMonitor(DA)
   D2spikes = SpikeMonitor(D2_L)
   D1_Lspikes = SpikeMonitor(D1_L)
   D1_NLspikes = SpikeMonitor(D1_NL)
   D1_WLspikes = SpikeMonitor(D1_WL) 
   
   win = [0.25, 0.25, 0.5]
   def calculate_FR(spike_monitor,integration_window,t): 
       bin_edges = np.arange((t-integration_window)/ms,t/ms,integration_window/4/ms)
       all_rates = []
       for i in range(0,len(bin_edges)-1):
           rate = np.size(np.where((spike_monitor.t > bin_edges[i]*ms)&(spike_monitor.t < bin_edges[i+1]*ms)))/n/integration_window
           all_rates.append(rate)
           smooth_rate = np.convolve(all_rates,win,mode='valid')
       return smooth_rate[len(smooth_rate)-1]*Hz
   

   def complex_DA_LTP(t,SpikeMon,SpikeMon2,SpikeMon3,SynapseMon,SynapseMon2):          
       DAind = np.where(SpikeMon2.t > (t - window)) 
       #SynapseMon.w = SynapseMon.w * 0.999
       if len(DAind[0]) > 5: # was DA released?
           CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
           PresynapticInd = []
           for i in range(0,len(np.unique(CortexIndex))):
               pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
               PresynapticInd.extend(pre[0])
           act_MSN = np.where(SpikeMon3 > MSN_High) # find MSNs that are in a "high" state
           high_synaptic_ind = np.concatenate([np.where(SynapseMon.j == x)[0] for x in act_MSN[0]]) # synaptic indicies projecting to high-state MSNs
           s2 = set(high_synaptic_ind) # set object for high_synpatic ind
           strengthen_synapse = [val for val in PresynapticInd if val in s2] # strengthen synapses that have MSNs in up state and cortical/DA input
           not_strengthen_synapse = list(set(PresynapticInd) - set(strengthen_synapse))# weaken synapses that have glutamate but not upstate
           SynapseMon.w[strengthen_synapse] +=   3*(SynapseMon.traceCon[strengthen_synapse] * mean(SynapseMon2.traceCon))     
           SynapseMon.w[not_strengthen_synapse] -= (SynapseMon.w[not_strengthen_synapse] - CortexD1_start)/np.abs(SynapseMon.w[not_strengthen_synapse]/CortexD1_start)
           SynapseMon.w = clip(SynapseMon.w, 0, wmax)
           
   def DA_LTP(t,window,SpikeMon,SpikeMon2,SpikeMon3,SynapseMon,SynapseMon2):          
       DAind = np.where(SpikeMon2.t > (t - window)) 
       #SynapseMon.w = SynapseMon.w * 0.999
       if len(DAind[0]) > 5: # was DA released?
           CortexIndex = SpikeMon.i[SpikeMon.t > defaultclock.t-window]  # index of cortical neurons that fired
           PresynapticInd = []
           for i in range(0,len(np.unique(CortexIndex))):
               pre = np.where(SynapseMon.i == np.unique(CortexIndex)[i]-1)
               PresynapticInd.extend(pre[0])
           act_MSN = np.where(SpikeMon3 > MSN_High) # find MSNs that are in a "high" state
           high_synaptic_ind = np.concatenate([np.where(SynapseMon.j == x)[0] for x in act_MSN[0]]) # synaptic indicies projecting to high-state MSNs
           s2 = set(high_synaptic_ind) # set object for high_synpatic ind
           strengthen_synapse = [val for val in PresynapticInd if val in s2] # strengthen synapses that have MSNs in up state and cortical/DA input
           not_strengthen_synapse = list(set(PresynapticInd) - set(strengthen_synapse))# weaken synapses that have glutamate but not upstate
           SynapseMon.w[strengthen_synapse] +=   1 
           SynapseMon.w[not_strengthen_synapse] -= (SynapseMon.w[not_strengthen_synapse] - CortexD1_start)/np.abs(SynapseMon.w[not_strengthen_synapse]/CortexD1_start)
           SynapseMon.w = clip(SynapseMon.w, 0, wmax)        

#network operations     
   rew_win=10*ms        
   @network_operation(dt=rew_win)
   def Reward(t):
       SNrL_rate = calculate_FR(SNrLspikes,integration_window,t)
       SNrNL_rate = calculate_FR(SNrNLspikes,integration_window,t)
       SNrWL_rate = calculate_FR(SNrWLspikes,integration_window,t)
       action = np.where(np.min([SNrL_rate,SNrNL_rate,SNrWL_rate]))
       if action[0] == 0:
           if SNrL_rate < SNr_thresh:
              DA.I += 2 # If a reward was recieved give DA
           else :
               DA.I = 0
       else :
            DA.I = 0          
           
   @network_operation(dt=rew_win)
   def calculate_LTP(t):
       DA_LTP(t,rew_win,CortexLspikes,DASpikes,D1_Lspikes,CortexL_D1L,DA_D1L)    

   @network_operation(dt=rew_win)
   def calculate_LTP2(t):
       DA_LTP(t,rew_win,CortexNLspikes,DASpikes,D1_NLspikes,CortexNL_D1NL,DA_D1NL)           
    
   @network_operation(dt=rew_win)    
   def calculate_LTP3(t):
       DA_LTP(t,rew_win,CortexWLspikes,DASpikes,D1_WLspikes,CortexWL_D1WL,DA_D1WL)    
  
    # Population monitors     
   CortexLpop = PopulationRateMonitor(CortexL)
   CortexNLpop = PopulationRateMonitor(CortexNL)
   CortexWLpop = PopulationRateMonitor(CortexWL)
   DApop = PopulationRateMonitor(DA)
   D1Lpop = PopulationRateMonitor(D1_L)
   D1NLpop = PopulationRateMonitor(D1_NL)
   D1WLpop = PopulationRateMonitor(D1_WL)
   GPePop = PopulationRateMonitor(GPe_L)
   SNrPop = PopulationRateMonitor(SNrL)
   SNrPopNL = PopulationRateMonitor(SNrNL)
   SNrPopWL = PopulationRateMonitor(SNrWL)
   STNpop = PopulationRateMonitor(STN)
   ThalamusPopL = PopulationRateMonitor(ThalamusL)
   ThalamusPopNL = PopulationRateMonitor(ThalamusNL)
   ThalamusPopWL = PopulationRateMonitor(ThalamusWL)  
    
   run(learn_duration,report='text',report_period=report_time)
    
   ##Compare population monitor to online smoothing 
   SNrL_binnedFR = [] 
   SNrNL_binnedFR = []  
   SNrWL_binnedFR = []       
   for i in np.arange(integration_window/ms,learn_duration/ms,rew_win/ms):
       SNrL_FR = calculate_FR(SNrLspikes,integration_window,i*ms)       
       SNrL_binnedFR.append(SNrL_FR)     
       SNrNL_FR = calculate_FR(SNrNLspikes,integration_window,i*ms)       
       SNrNL_binnedFR.append(SNrNL_FR)     
       SNrWL_FR = calculate_FR(SNrWLspikes,integration_window,i*ms)       
       SNrWL_binnedFR.append(SNrWL_FR)     
       
   figure()
   plot(SNrL_binnedFR,'r')
   plot(SNrNL_binnedFR,'b')
   plot(SNrWL_binnedFR,'g')
   plot(np.ones(len(SNrL_binnedFR))*SNr_thresh,'k--')
   legend('R2U')
   title('Online Calculation of SNr Firing rates')    
   plt.savefig(save_root + 'learnAction_SnrOnlinefiringRate.png')
   
   figure()
   plot(SNrPop.t/ms,SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(SNrPopNL.t/ms,SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(SNrPopWL.t/ms,SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('SNr Firing Rates')
   legend('R2U')
   plt.savefig(save_root + 'learnAction_SNrfiringRate.png')
   
   print 'Learned Action Results' 
   print np.str(np.sum(np.less(SNrL_binnedFR, SNr_thresh))) + '...rewarded action'
   print np.str(np.sum(np.less(SNrNL_binnedFR, SNr_thresh))) + '...unrewared action'
   print np.str(np.sum(np.less(SNrWL_binnedFR, SNr_thresh))) + '...unrewared action'
     
   avg_weights = [np.mean(CortexL_D1L.w), np.mean(CortexNL_D1NL.w), np.mean(CortexWL_D1WL.w)]
   print avg_weights
   
   if np.where(avg_weights == np.max(avg_weights))[0] == 0:
      print 'Learn action test passed'
   else:
          print 'Learn action test FAILED!!'
   
   figure()
   plot(DApop.t/ms,DApop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('DA Firing Rates')
   plt.savefig(save_root + 'learnAction_DAfiringRate.png')
    
   figure()
   plot(CortexLpop.t/ms,CortexLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(CortexNLpop.t/ms,CortexNLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(CortexWLpop.t/ms,CortexWLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('Cortex Firing Rates')
   legend('R2U')
   plt.savefig(save_root + 'learnAction_CortexfiringRate.png')
  
   figure()
   plot(D1Lpop.t/ms,D1Lpop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(D1NLpop.t/ms,D1NLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(D1WLpop.t/ms,D1WLpop.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('D1 Firing Rates')
   legend('R2U')   
   plt.savefig(save_root + 'learnAction_D1firingRate.png')
   

def test_synfire():
    Mgp = SpikeMonitor(CortexL) 
    CortexL.I = 0
    run(synfire_duration,report='text')
   
    figure()
    plot(Mgp.t/ms, 1.0*Mgp.i/group_size, '.')
    plot([0, synfire_duration/ms], np.arange(n_groups).repeat(2).reshape(-1, 2).T, 'k-')
    ylabel('group number')
    yticks(np.arange(n_groups))
    xlabel('time (ms)')
    title('Cortical Synfire for LevPress Channel')
    plt.savefig(save_root + 'SynfireTest.png')
   
   
def test_DA():
   CortexPop = PopulationRateMonitor(CortexL)
   DApop = PopulationRateMonitor(DA)
   D1pop = PopulationRateMonitor(D1_L)
   GPePop = PopulationRateMonitor(GPe_L)
   SNrPop = PopulationRateMonitor(SNrL)
   SNrPopNL = PopulationRateMonitor(SNrNL)
   SNrPopWL = PopulationRateMonitor(SNrWL)
   STNpop = PopulationRateMonitor(STN)
   ThalamusPopL = PopulationRateMonitor(ThalamusL)
   ThalamusPopNL = PopulationRateMonitor(ThalamusNL)
   ThalamusPopWL = PopulationRateMonitor(ThalamusWL)  

   DAvolts = StateMonitor(DA,('v'),record=True) 
   
   ta_DA = TimedArray([5,0,5],dt=DA_duration)
  
   @network_operation(dt=DA_duration)
   def update_sequence(t):
       DA.I = ta_DA(t)
       
   run(DA_duration*3,report='text',report_period=report_time)

   figure()
   plot(DAvolts.t,DAvolts.v[0])    

   figure()
   plot(SNrPop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   plot(SNrPopNL.smooth_rate(window='gaussian',width=binSize)/Hz,'b')
   plot(SNrPopWL.smooth_rate(window='gaussian',width=binSize)/Hz,'g')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('SNr Firing Rates')
   legend('R2U')
   plt.savefig(save_root + 'DAtest_SNrfiringRate.png')

   print  'DA firing rate=' + np.str(np.mean(DApop.rate))
   if np.mean(DApop.rate) > 0:
      print 'Test passed... DA neurons fire' 

   figure()
   plot(DApop.smooth_rate(window='gaussian',width=binSize)/Hz,'r')
   xlabel('Time(ms)')
   ylabel('Firing Rate')
   title('DA Firing Rates')
   plt.savefig(save_root + 'DAtest_DAfiringRate.png')