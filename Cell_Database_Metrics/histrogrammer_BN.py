#Created on June 26, 2017
#Modified on Aug 1, 2017
#Author: Anton Lunyov
#Pulls out specified data from 4D tensors generated for analyzing the effects of a,b,d,kLow on SFA, PIR and rheobase current

from pylab import *


##obtained from bubble_rheo_tensor.py, but the same for all (for consistency)
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


#Load all relevant data
transitionCurrents = load('TC_across_4.npy')
adaptation = load('adaptation_across_4.npy')
rheo = load('rheo_across_4.npy')



#H = 0.5, M = 0.3, L = 0.1
SFA_target = 0.46
rheo_target = 4.0
TC_target = -5.0

SFA_wind_arr = (0.1,0.45)
rheo_wind_arr = (1.0,3.0)
tc_wind_arr = (1.0,5.0)

SFA_names = ('N','B')
rheo_names = ('N','B')
TC_names = ('N','B')

for i in range(2):
    for j in range(2):
        for k in range(2):

            SFA_wind = SFA_wind_arr[i]
            rheo_wind = rheo_wind_arr[j]
            tc_wind = tc_wind_arr[k]
            
            SFA_up = SFA_target + SFA_wind
            SFA_low = SFA_target - SFA_wind

            rheo_up = rheo_target + rheo_wind
            rheo_low = rheo_target - rheo_wind

            TC_up = TC_target + tc_wind
            TC_low = TC_target - tc_wind

            inRange = (transitionCurrents < TC_up)*(transitionCurrents > TC_low)
            inRange = inRange * (adaptation < SFA_up)*(adaptation > SFA_low)
            inRange = inRange * (rheo < rheo_up)*(rheo > rheo_low)

            #Get the parameter values within
            aValsWithin = aVals[nonzero(inRange)[0]]
            bValsWithin = bVals[nonzero(inRange)[1]]
            dValsWithin = dVals[nonzero(inRange)[2]]
            kValsWithin = kVals[nonzero(inRange)[3]]

            #suptitle(str(float(len(aValsWithin))*100/10000) + "% of the models were explored out of all")
            suptitle(SFA_names[i] + rheo_names[j]+ str(TC_names[k]) + ": " + str(len(aValsWithin)) + " models were within range")
            subplot(221)
            hist(aValsWithin, bins = arange(aLow,aUp+aRes,aRes))
            xlabel('a')
            ylabel('Frequency of parameter')

            subplot(222)
            hist(bValsWithin, bins = arange(bLow,bUp+bRes,bRes))
            xlabel('b')
            ylabel('Frequency of parameter')

            subplot(223)
            hist(dValsWithin, bins = arange(dLow,dUp+dRes,dRes))
            xlabel('d')
            ylabel('Frequency of parameter')

            subplot(224)
            hist(kValsWithin, bins = arange(kLowLower,kLowUp+kRes,kRes))
            xlabel('k')
            ylabel('Frequency of parameter')

            show()
