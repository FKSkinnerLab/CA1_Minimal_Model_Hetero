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
SFA_target_arr = array([0.5,0.3,0.1])
SFA_names = ('H','M','L')
#H = 6.0 M = 4.0 L = 2.0
rheo_target_arr = array([6.0,4.0,2.0])
rheo_names = ('H','M','L')
#H = -10 M = -7 L = -4
TC_target_arr = array([-10.0,-7.0,-4.0])
TC_names = ('H','M','L')

SFA_wind = 0.1
rheo_wind = 1.0
tc_wind = 1.0

for i in range(3):
    for j in range(3):
        for k in range(3):
            
            SFA_target = SFA_target_arr[i]
            rheo_target = rheo_target_arr[j]
            TC_target = TC_target_arr[k]
            
            
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
