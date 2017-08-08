

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 10:05:50 2016

@author: elijah
"""

###############
# Basal Ganglia for action selection
# EP 170316
###############


#%% import and setup 
import multiprocessing 

import sys
sys.path.insert(0, 'C:\Users\elijah\Documents\GitHub\ModelingSNN\Brian2\BasalGanglia')
from basalGangliaTests import pop_firing, sequence, cortex_D1_action, test_DA, learn_action, MSN_SNr_inhibition

#learn_action()
#pop_firing()
#cortex_D1_action()
#MSN_SNr_inhibition()
#sequence()
#test_DA()

#%% Run and analyze
   
# will only work when run from the command line 
if __name__ == '__main__':
   p0 = multiprocessing.Process(name='p0',target = pop_firing, args = ())
   p1 = multiprocessing.Process(name='p1',target = sequence, args = ())
   p2 = multiprocessing.Process(name='p2',target = cortex_D1_action, args = ())
   p3 = multiprocessing.Process(name='p3',target = test_DA, args = ())
   p4 = multiprocessing.Process(name='p4',target = learn_action, args = ())
   p5 = multiprocessing.Process(name='p5',target = learn_action, args = ())
   
   p0.start()
   p1.start()
   p2.start()
   p3.start()   
   p4.start()
   p5.start()
   