from brian2 import *

start_scope()

tau = 20*ms
n = 100
b = 1.2 # constant current mean, the modulation varies
freq = 10*Hz
duration=1000*ms

eqs = '''
dv/dt = (-v + a * sin(2 * pi * freq * t) + b) / tau : 1
a : 1
'''

#reset='''
#v=0*mV
#vt += 3*mV
#'''

neurons = NeuronGroup(n, model=eqs, threshold='v > b+0.3', reset='v = 0.8',
                      method='euler')
neurons.v = 'rand()'
neurons.a = '0.05 + 0.7*i/n'
S = SpikeMonitor(neurons)
trace = StateMonitor(neurons, 'v', record=50)

run(duration)
subplot(211)
plot(S.t/ms, S.i, '.k')
xlabel('Time (ms)')
ylabel('Neuron index')
subplot(212)
plot(trace.t/ms, trace.v.T)
xlabel('Time (ms)')
ylabel('v')
show()


figure(2)
plot(neurons.v, S.count/duration)