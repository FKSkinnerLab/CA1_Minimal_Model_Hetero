# -*- coding: utf-8 -*-
#Created on June 2, 2017
#Modified on July 6, 2017
#Author: Anton Lunyov
#Generates a 4D tensor of SFA data as a function of a,b,d and kLow in excitatory cells based on the strongly adapting model

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

#Single neuron simulation - 50 neurons with varying current steps (for faster processing)
N = 50

#Contains ISI data
freq_initial = zeros(N)
freq_final = zeros(N)

#The values of the default current input - max and min. Sequence of events: Down-up-down
stepDown = 0.0 * pA
stepUp = 2.0 * pA

#Time breakdowns - the interval of applied current is between beforeStep and afterStep
beforeStep = 0.5 * second
afterStep = 1.5 * second

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

#Tensors to store all data obtained - 10x10x10x10
adaptation = zeros((len(aVals),len(bVals),len(dVals),len(kVals)))
adaptationRatio = zeros((len(aVals),len(bVals),len(dVals),len(kVals)))
#rheo = zeros((len(aVals),len(bVals),len(dVals),len(kVals)))

#Total iteration number - 0 to 10,000
count = 0

#Iterate over all parameters
for y in range(len(aVals)):
    for jj in range(len(bVals)):
        for ii in range(len(dVals)):
            for z in range(len(kVals)):
                
                #Print every 10th completed step
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
                Iapp = (t < beforeStep)*stepDown + (t >= beforeStep)*stepUp * i  + (t >= afterStep)*(-stepUp * i): amp

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
                #Cells_V = StateMonitor(CELLS,'v',record=True) #Don't record to speed up runtime
                Spiketimes = SpikeMonitor(CELLS)
				
		#Integration timecuts
                defaultclock.dt = 0.1*ms

                #Set duration and run simulation
                duration = 2 * second
                run(duration)

                for h in range(N):
                    #Retrieve the number of spikes for every neuron (and thus every current step)
                    spikes = Spiketimes.t[Spiketimes.i == h]
                    spikeNum = len(spikes)

                    #Record ISI based on conditions
                    if (spikeNum == 0):
                        freq_initial[h] = 0
                        freq_final[h] = 0
                    elif (spikeNum == 1):
                        freq_initial[h] = 1
                        freq_final[h] = 1
                    else:
                        #First and last ISI
                        freq_initial[h] = 1/(spikes[1] - spikes[0])
                        freq_final[h] = 1/(spikes[spikeNum - 1] - spikes[spikeNum - 2])

                #Re-create the implicit x array of currents to fit to
                x = range(0,int(stepUp*N/pA),int(stepUp/pA))
                #Fit initial and final adaptation data
                fitInitial = polyfit(x,freq_initial,1)
                slopeInitial = fitInitial[0]
                rheoInitial = -fitInitial[1]/slopeInitial

                fitFinal = polyfit(x,freq_final,1)
                slopeFinal = fitFinal[0]

                #Calculate relevant params
                SFA = slopeInitial - slopeFinal
                SFA_rat = (slopeInitial/slopeFinal)

                adaptation[y,jj,ii,z] = SFA
                adaptationRatio[y,jj,ii,z] = SFA_rat
		#Estimate of rheobase current - not implemented
                #rheo[y,jj,ii,z] = rheoInitial

#Save all data as numpy tensors
save('adaptation_across_4', adaptation)
save('adaptation_across_4_ratio',adaptationRatio)
print("Saved!")



