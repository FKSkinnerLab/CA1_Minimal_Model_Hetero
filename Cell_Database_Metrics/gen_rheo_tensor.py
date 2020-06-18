# -*- coding: utf-8 -*-
#Created on June 2, 2017
#Modified on July 6, 2017
#Author: Anton Lunyov
#This script creates a measure of rheobase current as a function of 4 parameters - a,b,d and kLow in excitatory cells based on the strongly adapting model

from brian2 import *
from pylab import *

#No compiling - remove if needed
prefs.codegen.target = "numpy"
#Don't include previous Brian objects
start_scope()

# Strongly adapting params (Izhikevitch)
Cm = 115.0 * pF # Membrane capacitance
vR = -61.8 * mV # Resting membrane potential
vPeak = 22.6 * mV # Spike cutoff
vT = -57.0 * mV # Instantaneous membrane potential
Ishift = 0.0 * pA # Rheobase current shift
Iother = 0.0 * pA # Defined at zero since input is not noisy

kLow = 0.1 * nS/mV
kHigh = 3.3  * nS/mV

a = 0.0012 /ms
d = 10.0 * pA
b = 3.0 * nS
c = -65.8 * mV

#Single neuron simulation - 100 neurons with varying current steps
N = 100

#The resolution of the current steps tested. In this case, the current goes from -25 to 25 by 0.5 increments
step = 0.5 * pA
offset = -(step * N)/2

#Time breakdowns - the interval of applied current is between beforeStep and afterStep
beforeStep = 0.5 * second

# a value iterations
aRes = 2.4E-4 # 1.0E-4
aUp = 2.4E-3
aLow = 0.0E-3
aVals = arange(aLow,aUp,aRes)

# klow iterations
kRes = 0.02 # 0.01
kLowUp = 0.2
kLowLower = 0.0
kVals = arange(kLowLower,kLowUp,kRes)

# b iters
bRes = 0.6 #0.2
bUp = 6.0
bLow = 0.0
bVals = arange(bLow,bUp,bRes)

# d iters
dRes = 2 # 0.5
dUp = 20
dLow = 0
dVals = arange(dLow,dUp,dRes)

#Tensor to store data (10x10x10x10)
rheo = zeros((len(aVals),len(bVals),len(dVals),len(kVals)))

#Total iteration number - 0 to 10,000
count = 0

#Iterate over all parameters
for y in range(len(aVals)):
    for jj in range(len(bVals)):
        for ii in range(len(dVals)):
            for z in range(len(kVals)):
			
                #Print every 100th completed step
                count = count + 1
                if (count % 100 == 0):
                    print(count)
                
		#Find parameter values for all neurons at current iteration
                a = aVals[y] /ms
                kLow = kVals[z] * nS/mV
                d = dVals[ii] * pA
                b = bVals[jj] * nS

                #Izhikevitch-Skinner cell equations
                cellEqs = """
                #Excitatory term
                Iapp = (t > beforeStep)*((step*(i-1)) + offset): amp

                k = (v<vT)*kLow+(v>=vT)*kHigh : siemens/volt
                du/dt = a*(b*(v-vR)-u) : amp
                dv/dt = (k*(v-vR)*(v-vT)+Iapp-u)/Cm : volt 
                """

                #Create cells
                CELLS = NeuronGroup(N, model=cellEqs, reset ="v = c; u += d" , threshold="v>=vPeak", method="euler")

                #set initial conditions
                #Resting potential used as base
                CELLS.v = vR


                #Record data
                #Cells_V = StateMonitor(CELLS,'v',record=True)
                Spiketimes = SpikeMonitor(CELLS)
				
		#Integration timecuts
                defaultclock.dt = 0.1*ms

                #Set duration and run simulation
                duration = 2 * second

                run(duration)

                #For every neuron
                for h in range(N):
                    #Retrieve the number of spikes
                    spikes = Spiketimes.t[Spiketimes.i == h]
                    spikes = spikes[spikes < 1.5*second]
                    spikeNum = len(spikes)

                    #Record rheobase current
                    if (spikeNum > 0):
                        current = (step*h) + offset
                        rheo[y,jj,ii,z] = current/pA
                        break
		    #If no spike occured (end is reached), record maximum current (serves as a cutoff, 25 and above are all 25)
                    elif (h == N - 1):
                        rheo[y,jj,ii,z] = (step*h) + offset

#Save numpy array as file for 4D_pullout.py
save('rheo_across_4',rheo)
print("Saved")



