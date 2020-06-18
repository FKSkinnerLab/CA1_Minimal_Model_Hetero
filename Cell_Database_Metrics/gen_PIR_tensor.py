# -*- coding: utf-8 -*-
#Created on June 2, 2017
#Modified on Aug 15, 2017
#Author: Anton Lunyov
#This script creates a measure of post-inhibitory rebound as a function of 4 parameters - a,b,d and kLow in excitatory cells based on the strongly adapting model

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

#Single neuron simulation - 50 neurons with varying current steps
N = 50

#The values of the default current input - max and min. Down-up-down
stepDown = 0.0 * pA
stepUp = -0.5 * pA

#Time breakdowns - the interval of applied current is between beforeStep and afterStep. Duration of 1 second
beforeStep = 0.5 * second
afterStep = 1.5 * second

# a value iterations - resolution is 10
aRes = 2.4E-4 # 1.0E-4
aUp = 2.4E-3
aLow = 0.0E-3
aVals = arange(aLow,aUp,aRes)

# klow iterations - resolution is 10
kRes = 0.02 # 0.01
kLowUp = 0.2
kLowLower = 0.0
kVals = arange(kLowLower,kLowUp,kRes)

# b iters - resolution is 10
bRes = 0.6 #0.2
bUp = 6.0
bLow = 0.0
bVals = arange(bLow,bUp,bRes)

# d iters - resolution is 10
dRes = 2 # 0.5
dUp = 20
dLow = 0
dVals = arange(dLow,dUp,dRes)

#Tensor to store data
transitionCurrents = zeros((len(aVals),len(bVals),len(dVals),len(kVals)))
#Total iteration number - 0 to 10,000
count = 0

#For each of the parameters
for y in range(len(aVals)):
    for z in range(len(kVals)):
        for ii in range(len(dVals)):
            for jj in range(len(bVals)): 

                #Keep track of iteration number and print out every xth
                count = count + 1
                if (count % 1000 == 0):
                    print(count)

                #Find parameter values at current iteration
                a = aVals[y] /ms
                kLow = kVals[z] * nS/mV
                d = dVals[ii] * pA
                b = bVals[jj] * nS

                #Izhikevitch-Skinner cell equations
                cellEqs = """
                #Excitatory term
				#The excitatory term to each cell depends on the cell number - that way, implemented more efficiently as current steps down
                Iapp = (t < beforeStep)*stepDown + (t >= beforeStep)*stepUp * i  + (t >= afterStep)*(-stepUp * i): amp

                k = (v<vT)*kLow+(v>=vT)*kHigh : siemens/volt
                du/dt = a*(b*(v-vR)-u) : amp
                dv/dt = (k*(v-vR)*(v-vT)+Iapp-u)/Cm : volt 
                """

                #Create cells based on models
                CELLS = NeuronGroup(N, model=cellEqs, reset ="v = c; u += d" , threshold="v>=vPeak", method="euler")

                #set initial conditions
                #Resting potential used as base
                CELLS.v = vR


                #Record plotted data
                Cells_V = StateMonitor(CELLS,'v',record=True)
                Spiketimes = SpikeMonitor(CELLS)
				
		#Integration timecuts 
                defaultclock.dt = 0.1*ms

                #Set duration and run simulation
                duration = 2 * second

                run(duration)

                #These neurons are spiking before the membrane potential is decreased, and are thus ineligible for PIR analysis (in my limited expertise)
                falseSpikes = len(Spiketimes.t[Spiketimes.t < 1.5*second]/second)
                
                #If was not spiking before,
                if (falseSpikes == 0):
                        spikes = zeros(N)
                        for h in range(N):
                                #Record the number of spikes after the hyperpolarization is over
                                hSpikes = Spiketimes.t[Spiketimes.i == h]
                                hSpikes = hSpikes[hSpikes >= 1.5*second]/second
                                spikes[h] = len(hSpikes)

                        #Check where the transition occurs, put one zero so array is not empty (checks are easier)
                        transitions = zeros(1)
                        
                        #If at zero point, save that as transition current and give warning
                        if (spikes[0] > 1):
                                transitions = append(transitions,0.0)
                                print("Transition current of 0 occured")
                        
			#Otherwise find the exact point where the transition occurs. Choose latter of the two as the point
                        for k in range(1,len(spikes)):
                                if ((int(spikes[k]) > 0) & (int(spikes[k - 1]) == 0)):
                                        transitions = append(transitions,k*stepUp/pA)
                        
                        #Eliminate the zero inserted in beginning, and record number of transitions
                        transitionNum = len(transitions) - 1
                        
                        #Invalid TC = +1 pA
                        #If more than one transitions occured, the effect is invalid - model breaks down
                        if (transitionNum > 1):
                                transitionCurrents[y,jj,ii,z] = +1
                        #No transition occured - current required is not enough, thus invalid
                        elif (transitionNum == 0):
                                transitionCurrents[y,jj,ii,z] = +1
                        elif (transitionNum == 1):
                                transitionCurrents[y,jj,ii,z] = transitions[1]
                                #print(transitions[1])
                #If cell firing before, invalid
                else:
                        transitionCurrents[y,jj,ii,z] = +1


#Save numpy array as file to be used in 4D_pullout.py
save('TC_across_4',transitionCurrents)
print("Saved!")
