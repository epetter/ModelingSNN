# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 10:49:20 2016

@author: Elijah
"""
import numpy as np

## BG analysis

rangez = 10000
# Extract
IntThal = ThalamusPop.rate[0:rangez]/n/Hz
IntThalT = ThalamusPop.t[0:rangez]/ms
EndThal = ThalamusPop.rate[len(ThalamusPop.rate)-rangez:len(ThalamusPop.rate)]/n/Hz
EndThalT = ThalamusPop.t[len(ThalamusPop.t)-rangez:len(ThalamusPop.t)]/ms

IntAction = ActionPop.rate[0:rangez]/n/Hz
IntActionT = ActionPop.t[0:rangez]/ms
EndAction = ActionPop.rate[len(ActionPop.rate)-rangez:len(ActionPop.rate)]/n/Hz
EndActionT = ActionPop.t[len(ActionPop.t)-rangez:len(ActionPop.t)]/ms

IntNoAct = NoActionPop.rate[0:rangez]/n/Hz
IntNoActT = NoActionPop.t[0:rangez]/ms
EndNoActT = NoActionPop.t[len(NoActionPop.t)-rangez:len(NoActionPop.t)]/ms
EndNoAct = NoActionPop.rate[len(NoActionPop.rate)-rangez:len(ActionPop.rate)]/n/Hz
 
IntWrongActT = WrongActionPop.t[0:rangez]/ms
IntWrongAct =  WrongActionPop.rate[0:rangez]/n/Hz
EndWrongActT = WrongActionPop.t[len(WrongActionPop.t)-rangez:len(WrongActionPop.t)]/ms
EndWrongAct = WrongActionPop.rate[len(WrongActionPop.rate)-rangez:len(WrongActionPop.rate)]/n/Hz

# Plot
figure(figsize=(10, 20))
subplot(421)
plot(IntThalT,IntThal)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Initial Thalamus Population Activity')
subplot(422)
plot(EndThalT,EndThal)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Ending Thalamus Population Activity')
subplot(423)
plot(IntActionT,IntAction)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Initial Action #1 Activity')
subplot(424)
plot(EndActionT,EndAction)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Ending Action #1 Activity')
subplot(425)
plot(IntNoActT,IntNoAct)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Initial Action #2 Activity')
subplot(426)
plot(EndNoActT,EndNoAct)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Ending Action #2 Activity')
subplot(427)
plot(IntWrongActT,IntWrongAct)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Initial Action #3 Activity')
subplot(428)
plot(EndWrongActT,EndWrongAct)
ylim(0,50)
xlabel('Time(ms)')
ylabel('Rate(Hz)')
title('Ending Action #3 Activity')

# Calcukate continuous correlation

ThalAct = [IntThal,IntAction]
EndThalAct = [EndThal,EndAction]
IntActionCorr = np.corrcoef(ThalAct)
EndActionCorr = np.corrcoef(EndThalAct)

ThalNoAct = [IntThal,IntNoAct]
EndThalNoAct = [EndThal,EndNoAct]
ThalNoActCorr = np.corrcoef(ThalNoAct)
EndThalNoActCorr = np.corrcoef(EndThalNoAct)

ThalWrongAct = [IntThal,IntWrongAct]
EndThalWrongAct = [EndThal,EndWrongAct]
ThalWrongActCorr = np.corrcoef(ThalWrongAct)
EndThalWrongActCorr = np.corrcoef(EndThalWrongAct)
