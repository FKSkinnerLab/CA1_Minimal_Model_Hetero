# -*- coding: utf-8 -*-
#Created on July 10, 2017
#Modified on August 1, 2017
#Author: Anton Lunyov

from brian2 import *
from pylab import *
from bursts_brian2 import find_bursts,find_distn,make_dist_plots,make_binwidth_plot
from determine_stats_brian2 import bursting_stats, make_stats_plots, make_stats_plots
from write_stats_brian2 import write_stats
from scipy.signal import find_peaks_cwt
import os

#No compiling to speed up process
prefs.codegen.target = "numpy"

#Set random seed for replication
seed(71)

#Input for tensor window sizes
arg_sfa_wind = sys.argv[1]
arg_rheo_wind = sys.argv[2]
arg_pir_wind = sys.argv[3]

arg_sfa_tar = sys.argv[4]
arg_rheo_tar = sys.argv[5]
arg_pir_tar = sys.argv[6]

#Connectivities
pyr_pv = 0.02
pv_pyr = 0.3
pyr_pyr = 0.01
pv_pv = 0.12

#Cell numbers
N_pv = 500
N_pyr = 10000

#Integration timestep
defaultclock.dt = 0.04*ms

######################Model parameters (intrinsic parameters in Katie's code)

#PYR cells' strongly adapting model params (Izhikevitch)
Cm_pyr = 115.0 * pF # Membrane capacitance
vR_pyr = -61.8 * mV # Resting membrane potential
vPeak_pyr = 22.6 * mV # Spike cutoff
vT_pyr = -57.0 * mV # Instantaneous membrane potential
Ishift_pyr = 0.0 * pA # Rheobase current shift
Iother_pyr = 0.0 * pA
kLow_pyr = 0.1 * nS/mV
kHigh_pyr = 3.3  * nS/mV
a_pyr = 0.0012 /ms
d_pyr = 10.0 * pA
b_pyr = 3.0 * nS
c_pyr = -65.8 * mV

#Obtained from PYR_PV_Network_compact.py
#PV cells' model parameters
Cm_pv = 90.0 * pF
vR_pv = -60.6 * mV 
vPeak_pv = -2.5 * mV
vT_pv = -43.1 * mV
Ishift_pv = 0.0 * pA
Iother_pv = 0.0 * pA
kLow_pv = 1.7 * nS/mV
kHigh_pv = 14.0  * nS/mV
a_pv = 0.1 /ms
d_pv = 0.1 * pA
b_pv = -0.1 * nS
c_pv = -67.0 * mV


######################Synaptic parameters

Erev_e = -15.0 * mV #Reversal potential of excitatory cell
Erev_i = -85.0 * mV #Reversal potential of inhibitory cell

#PYR to PYR synapses
g_pyr_pyr = 0.094 * nS # biological 0.008-5 nS
alpha_pyr_pyr = 2000.0 /second
beta_pyr_pyr = 333.33 /second
s_inf_pyr_pyr = alpha_pyr_pyr/(alpha_pyr_pyr+beta_pyr_pyr)
tau_s_pyr_pyr = 1.0/(alpha_pyr_pyr+beta_pyr_pyr)

#PV to PYR synapses
g_pv_pyr = 8.7 * nS
alpha_pv_pyr = 3333.0/second
beta_pv_pyr = 286.0/second
s_inf_pv_pyr = alpha_pv_pyr/(alpha_pv_pyr+beta_pv_pyr)
tau_s_pv_pyr = 1.0/(alpha_pv_pyr+beta_pv_pyr)

#PYR to PV synapses
g_pyr_pv = 3.0 * nS
alpha_pyr_pv = 2710.0 / second
beta_pyr_pv = 483.0 / second
s_inf_pyr_pv = alpha_pyr_pv/(alpha_pyr_pv+beta_pyr_pv)
tau_s_pyr_pv = 1.0/(alpha_pyr_pv+beta_pyr_pv)

#PV to PV synapses
g_pv_pv = 3.0 * nS
alpha_pv_pv = 3700.0 / second
beta_pv_pv = 560.0 / second
s_inf_pv_pv = alpha_pv_pv/(alpha_pv_pv+beta_pv_pv)
tau_s_pv_pv = 1.0/(alpha_pv_pv+beta_pv_pv)

#OU parameters
tau_excit=2.73*ms
ge_mean=0.0
ge_SD=0.6

#############################################################################
#Set parameters on PYR cells one by one to make the population heterogenous

#Load all relevant data
transitionCurrents = load('TC_across_4.npy')
adaptation = load('adaptation_across_4.npy')
rheo = load('rheo_across_4.npy')

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

#Values obtained at default parameter values
SFA_target = float(arg_sfa_tar)#0.46 
rheo_target = float(arg_rheo_tar)#4.0
TC_target = float(arg_pir_tar)#-5.0

SFA_window = float(arg_sfa_wind)#0.45
rheo_window = float(arg_rheo_wind)#3.0
TC_window = float(arg_pir_wind)#5.0

SFA_up = SFA_target + SFA_window
SFA_low = SFA_target - SFA_window

rheo_up = rheo_target + rheo_window
rheo_low = rheo_target - rheo_window

TC_up = TC_target + TC_window
TC_low = TC_target - TC_window

inRange = (transitionCurrents < TC_up)*(transitionCurrents > TC_low)
inRange = inRange * (adaptation < SFA_up)*(adaptation > SFA_low)
inRange = inRange * (rheo < rheo_up)*(rheo > rheo_low)

#Get the parameter values within
aValsWithin = aVals[nonzero(inRange)[0]]
bValsWithin = bVals[nonzero(inRange)[1]]
dValsWithin = dVals[nonzero(inRange)[2]]
kValsWithin = kVals[nonzero(inRange)[3]]

#Length of the allowed parameter slice
lenOfAPS = len(aValsWithin)
print(str(lenOfAPS) + " models picked out of 10,000")

aSet = zeros(N_pyr)
bSet = zeros(N_pyr)
dSet = zeros(N_pyr)
kSet = zeros(N_pyr)

h = array(rand(N_pyr)*N_pyr,dtype = int)

for c in range(N_pyr):
    #Index of model parameters
    #h = #rand(0,lenOfAPS)#c % lenOfAPS
    q = h[c] % lenOfAPS
    #print(b)
    aSet[c] = aValsWithin[q]
    bSet[c] = bValsWithin[q]
    dSet[c] = dValsWithin[q]
    kSet[c] = kValsWithin[q]


#Izhikevitch-Skinner cell equations for E and I cells
cellEqs_pyr = '''
#Model params explicitly defined to allow them to be modified
a : second**-1
b : siemens
d : amp
kLow : siemens/volt

#Synaptic input
Isyn_pv_pyr = g_pv_pyr*(v-Erev_i)*s_sum_pv_pyr : amp
s_sum_pv_pyr : 1

Isyn_pyr_pyr = g_pyr_pyr*(v-Erev_e)*s_sum_pyr_pyr : amp
s_sum_pyr_pyr : 1

k = (v<vT_pyr)*kLow+(v>=vT_pyr)*kHigh_pyr : siemens/volt
du/dt = a*(b*(v-vR_pyr)-u) : amp
dv/dt = ((k*(v-vR_pyr)*(v-vT_pyr)) + Iext -Isyn_pv_pyr - Isyn_pyr_pyr -u)/Cm_pyr : volt

#Excitatory term - OU
Iext = (gext*nS)*(v-Erev_e) : amp
dgext/dt = -(gext-ge_mean)/tau_excit + sqrt(2*(ge_SD)**2/tau_excit)*xi : 1
'''
#cellEqs_pyr+=OrnsteinUhlenbeck('genoise',mu=ge_mean*nS,sigma=ge_SD*nS,tau=tau_excit)
#Iext=genoise*(v-Ee) : amp

cellEqs_pv = '''
#Excitatory term
Iext : amp

#Synaptic input
Isyn_pyr_pv = g_pyr_pv*(v-Erev_e)*s_sum_pyr_pv : amp
s_sum_pyr_pv : 1

Isyn_pv_pv = g_pv_pv*(v-Erev_i)*s_sum_pv_pv : amp
s_sum_pv_pv : 1

k = (v<vT_pv)*kLow_pv+(v>=vT_pv)*kHigh_pv : siemens/volt
du/dt = a_pv*(b_pv*(v-vR_pv)-u) : amp
dv/dt = ((k*(v-vR_pv)*(v-vT_pv)) + Iext - Isyn_pyr_pv - Isyn_pv_pv -u)/Cm_pv : volt 
'''

#Create cells based on models above
PYR = NeuronGroup(N_pyr, model=cellEqs_pyr, reset ="v = c_pyr; u += d" , threshold="v>=vPeak_pyr", method="euler")
PV = NeuronGroup(N_pv, model=cellEqs_pv, reset ="v = c_pv; u += d_pv" , threshold="v>=vPeak_pv", method="euler")
    

#################################################################################
#############Synaptic models

#PYR-PYR
synModel_pyr_pyr = '''
ds0/dt=-(s0-s_inf_pyr_pyr)/(tau_s_pyr_pyr) :1 
ds1/dt=-s1*beta_pyr_pyr :1 
s_tot=clip(s0,0,s1) :1
s_sum_pyr_pyr_post=s_tot:1(summed)
'''
synPre_pyr_pyr = '''
s0=s1;
s1=(s_inf_pyr_pyr*(1-exp(-0.001*second/tau_s_pyr_pyr))+s1*exp(-0.001*second/tau_s_pyr_pyr))*exp(0.001*beta_pyr_pyr*second)'''

SYN_PYR_PYR=Synapses(PYR,PYR,model=synModel_pyr_pyr, pre=synPre_pyr_pyr, method='rk2')

#PYR-PV
synModel_pyr_pv = '''
ds0/dt=-(s0-s_inf_pyr_pv)/(tau_s_pyr_pv) :1 
ds1/dt=-s1*beta_pyr_pv :1 
s_tot=clip(s0,0,s1) :1
s_sum_pyr_pv_post=s_tot:1(summed)
'''
synPre_pyr_pv = '''
s0=s1;
s1=(s_inf_pyr_pv*(1-exp(-0.001*second/tau_s_pyr_pv))+s1*exp(-0.001*second/tau_s_pyr_pv))*exp(0.001*beta_pyr_pv*second)'''

SYN_PYR_PV=Synapses(PYR,PV,model=synModel_pyr_pv,pre=synPre_pyr_pv, method='rk2')


#PV-PYR
synModel_pv_pyr = '''
ds0/dt=-(s0-s_inf_pv_pyr)/(tau_s_pv_pyr) :1 
ds1/dt=-s1*beta_pv_pyr :1 
s_tot=clip(s0,0,s1) :1
s_sum_pv_pyr_post=s_tot:1(summed)
'''
synPre_pv_pyr = '''
s0=s1;
s1=(s_inf_pv_pyr*(1-exp(-0.001*second/tau_s_pv_pyr))+s1*exp(-0.001*second/tau_s_pv_pyr))*exp(0.001*beta_pv_pyr*second)'''

SYN_PV_PYR=Synapses(PV,PYR,model=synModel_pv_pyr,pre=synPre_pv_pyr, method='rk2')

#PV-PV
synModel_pv_pv = '''
ds0/dt=-(s0-s_inf_pv_pv)/(tau_s_pv_pv) :1 
ds1/dt=-s1*beta_pv_pv :1 
s_tot=clip(s0,0,s1) :1
s_sum_pv_pv_post=s_tot:1(summed)
'''
synPre_pv_pv = '''
s0=s1;
s1=(s_inf_pv_pv*(1-exp(-0.001*second/tau_s_pv_pv))+s1*exp(-0.001*second/tau_s_pv_pv))*exp(0.001*beta_pv_pv*second)'''

SYN_PV_PV=Synapses(PV,PV,model=synModel_pv_pv, pre=synPre_pv_pv, method='rk2')


#set initial conditions
#-55 to -65 mV used as base

PYR.v = ((rand(N_pyr))*10.0 - 65.0)*mV
PYR.s_sum_pv_pyr = 0.0 
PYR.s_sum_pyr_pyr = 0.0

PV.v = ((rand(N_pv))*10.0 -  65.0)*mV
PV.s_sum_pyr_pv = 0.0 
PV.s_sum_pv_pv = 0.0


#Connect cells with synpases with a probability p
SYN_PYR_PYR.connect('i!=j', p=pyr_pyr)
SYN_PV_PYR.connect('i!=j',p=pv_pyr)
SYN_PYR_PV.connect('i!=j',p=pyr_pv)
SYN_PV_PV.connect('i!=j',p=pv_pv)

initial_values = {'a': aSet * 1/ms ,
                  'b': bSet * nS,
                  'd': dSet * pA,
                  'kLow': kSet *nS/mV}
#Make homogenous
PYR.set_states(initial_values)


#Record plotted data
N_pyr_mon=100
N_pv_mon=50

PYR_mon = PYR[0:N_pyr_mon]
PV_mon = PV[0:N_pv_mon]

Cells_V = StateMonitor(PYR_mon,'v',record=True)
Spiketimes_PYR = SpikeMonitor(PYR, record=True)
Spiketimes_PV = SpikeMonitor(PV, record=True)

PV_Esyn = StateMonitor(PV_mon,'Isyn_pyr_pv',record=True)
PV_Isyn = StateMonitor(PV_mon,'Isyn_pv_pv',record=True)
PYR_Esyn = StateMonitor(PYR_mon,'Isyn_pyr_pyr',record=True)
PYR_Isyn = StateMonitor(PYR_mon,'Isyn_pv_pyr',record=True)


#Set duration and run simulation
duration_num = 10
duration = duration_num * second 
run(duration)

#Analysis starts here
print("Simulation complete")
print("Starting analysis")

transient = 5000

#Array of average voltages by timestep
AvgVoltage = sum(Cells_V.v, axis = 0)

##calculate fourier transform
##Code adapted from PYR_bursting_networkx.py
###Ignore transient stage
skip_time=int((transient*ms)/defaultclock.dt)
network_fft = abs(fft(AvgVoltage[skip_time:])/len(AvgVoltage[skip_time:]))
length = len(network_fft)
freq = fftfreq(length, d=defaultclock.dt)
 
l = array(network_fft[1:length/2])
max_index = l.argmax()
max_freq=freq[max_index+1]
max_pwr=network_fft[max_index+1]

print("Max pwr: " + str(max_pwr) + " max freq: " +str(max_freq))
#print(Spiketimes_PYR.i)



#Adapted from PYR_PV_network_compact.py
[bin_s,bin_ms,binpts,totspkhist,totspkdist_smooth,dist_thresh,totspkhist_list,
         thresh_plot,binpt_ind,lgbinwidth,numlgbins,
         intraburst_bins,interburst_bins,intraburst_bins_ms,interburst_bins_ms,intraburst_time_ms,interburst_time_ms,
         num_intraburst_bins,num_interburst_bins,bin_ind_no_trans,intrabin_ind_no_trans,
         interbin_ind_no_trans] = find_bursts(duration, defaultclock.dt, transient, N_pyr,Spiketimes_PYR.t,Spiketimes_PYR.i,max_freq/Hz)


path = os.path.dirname(os.path.realpath(__file__))
foldername = path + '/Output'
#SFA, Rheo, PIR
filename = "[" + str(SFA_target - SFA_window) + "_" + str(SFA_target + SFA_window)
filename = filename + ',' + str(rheo_target - rheo_window) + '_' + str(rheo_target + rheo_window)
filename = filename + ',' + str(TC_target - TC_window) + '_' + str(TC_target + TC_window) + ']_dur' + str(duration_num)


g_pyrpyr_str = str(g_pyr_pyr)
ge_mean_str = str(ge_mean)
ge_SD_str = str(ge_SD)
celltype = ('PYR','PV')
num_popln = 2
M_t_list = (Spiketimes_PYR.t,Spiketimes_PV.t)
M_i_list = (Spiketimes_PYR.i,Spiketimes_PV.i)
N_list=[N_pyr,N_pv]


#plot distns and binwidths, and return binwidth_ms,ctr_of_bin_ms, and updated bin_ms & bin_s
#Order matters
make_dist_plots(duration,transient,num_popln,lgbinwidth,bin_s,bin_ms,binpts,binpt_ind,totspkdist_smooth,thresh_plot,
                dist_thresh,intraburst_time_ms,interburst_time_ms,max_freq,max_pwr,g_pyrpyr_str,ge_mean_str,
                ge_SD_str,celltype,filename,foldername)

[binwidth_ms,intrabinwidth_ms,interbinwidth_ms,ctr_of_bin_ms,ctr_of_intrabin_ms,
ctr_of_interbin_ms,numbins,bin_ms,bin_s]=make_binwidth_plot(duration,bin_s,bin_ms,intraburst_bins_ms,interburst_bins_ms,
                                                            transient,filename,foldername)

##make_dist_plots(duration,num_popln,lgbinwidth,bin_s,bin_ms,binpts,binpt_ind,totspkdist_smooth,thresh_plot,
##                dist_thresh,intraburst_time_ms,interburst_time_ms,max_freq,max_pwr,g_pyrpyr_str,ge_mean_str,
##                ge_SD_str,celltype,filename,foldername)
##
##[binwidth_ms,intrabinwidth_ms,interbinwidth_ms,ctr_of_bin_ms,ctr_of_intrabin_ms,
##ctr_of_interbin_ms,numbins,bin_ms,bin_s]=make_binwidth_plot(duration,bin_s,bin_ms,intraburst_bins_ms,interburst_bins_ms,
##                                                            transient,filename,foldername)

###Save
##plot(totspkdist_smooth,totspkhist,'.')
##xlabel("x")
##ylabel("Histogram")
##savefig(foldername+"/"+filename+"_PYR_hist.png")
##close()

print("Started stats")

#initialize lists
totspkhist=[[] for i in xrange(num_popln)]  # or can do: [[]]*num_popln   
totspkdist_smooth=[[] for i in xrange(num_popln)]
totspkhist_list=[[] for i in xrange(num_popln)]
indices=[[] for i in xrange(num_popln)] 
num_cells_w_FBF=[[] for i in xrange(num_popln)]
fullburst_freq_nonzero=[[] for i in xrange(num_popln)]
fullburst_freq_nonzero_no_trans=[[] for i in xrange(num_popln)]
#bin_ind_no_trans=[[] for i in xrange(num_popln)] 
num_cells_spk_per_bin=[[] for i in xrange(num_popln)]
num_cells_spk_per_bin_no_trans=[[] for i in xrange(num_popln)]
avg_num_cells_spk=[[] for i in xrange(num_popln)]
std_num_cells_spk=[[] for i in xrange(num_popln)]
num_spks_per_cell_per_bin=[[] for i in xrange(num_popln)]
num_spks_per_cell_per_bin_no_trans=[[] for i in xrange(num_popln)]
avg_num_spks_per_cell=[[] for i in xrange(num_popln)]
std_num_spks_per_cell=[[] for i in xrange(num_popln)]
tot_num_spks_per_bin=[[] for i in xrange(num_popln)]
tot_num_spks_per_bin_no_trans=[[] for i in xrange(num_popln)]
avg_tot_num_spks=[[] for i in xrange(num_popln)]
std_tot_num_spks=[[] for i in xrange(num_popln)]
avg_fullburst_freq_per_bin_tot=[[] for i in xrange(num_popln)]
std_fullburst_freq_per_bin_tot=[[] for i in xrange(num_popln)]
avg_fullburst_freq_per_cell_tot=[[] for i in xrange(num_popln)]
std_fullburst_freq_per_cell_tot=[[] for i in xrange(num_popln)]
avg_fullburst_freq_per_bin=[[] for i in xrange(num_popln)]
std_fullburst_freq_per_bin=[[] for i in xrange(num_popln)]
avg_FBF_per_bin_list=[[] for i in xrange(num_popln)]
std_FBF_per_bin_list=[[] for i in xrange(num_popln)]
avg_FBF_per_cell_list=[[] for i in xrange(num_popln)]
std_FBF_per_cell_list=[[] for i in xrange(num_popln)]
intraindices=[[] for i in xrange(num_popln)] 
num_cells_w_intraBF=[[] for i in xrange(num_popln)]
intraburst_freq_nonzero=[[] for i in xrange(num_popln)]
intraburst_freq_nonzero_no_trans=[[] for i in xrange(num_popln)]
#intrabin_ind_no_trans=[[] for i in xrange(num_popln)] 
num_cells_spk_per_intrabin=[[] for i in xrange(num_popln)]
num_cells_spk_per_intrabin_no_trans=[[] for i in xrange(num_popln)]
avg_num_cells_spk_intra=[[] for i in xrange(num_popln)]
std_num_cells_spk_intra=[[] for i in xrange(num_popln)]
num_spks_per_cell_per_intrabin=[[] for i in xrange(num_popln)]
num_spks_per_cell_per_intrabin_no_trans=[[] for i in xrange(num_popln)]
avg_num_spks_per_cell_intra=[[] for i in xrange(num_popln)]
std_num_spks_per_cell_intra=[[] for i in xrange(num_popln)]
tot_num_spks_per_intrabin=[[] for i in xrange(num_popln)]
tot_num_spks_per_intrabin_no_trans=[[] for i in xrange(num_popln)]
avg_tot_num_spks_intra=[[] for i in xrange(num_popln)]
std_tot_num_spks_intra=[[] for i in xrange(num_popln)]
avg_intraburst_freq_per_bin_tot=[[] for i in xrange(num_popln)]
std_intraburst_freq_per_bin_tot=[[] for i in xrange(num_popln)]
avg_intraburst_freq_per_cell_tot=[[] for i in xrange(num_popln)]
std_intraburst_freq_per_cell_tot=[[] for i in xrange(num_popln)]
avg_intraburst_freq_per_bin=[[] for i in xrange(num_popln)]
std_intraburst_freq_per_bin=[[] for i in xrange(num_popln)]
avg_intraBF_per_bin_list=[[] for i in xrange(num_popln)]
std_intraBF_per_bin_list=[[] for i in xrange(num_popln)]
avg_intraBF_per_cell_list=[[] for i in xrange(num_popln)]
std_intraBF_per_cell_list=[[] for i in xrange(num_popln)]
interindices=[[] for i in xrange(num_popln)] 
num_cells_w_interBF=[[] for i in xrange(num_popln)]
interburst_freq_nonzero=[[] for i in xrange(num_popln)]
interburst_freq_nonzero_no_trans=[[] for i in xrange(num_popln)]
#interbin_ind_no_trans=[[] for i in xrange(num_popln)] 
num_cells_spk_per_interbin=[[] for i in xrange(num_popln)]
num_cells_spk_per_interbin_no_trans=[[] for i in xrange(num_popln)]
avg_num_cells_spk_inter=[[] for i in xrange(num_popln)]
std_num_cells_spk_inter=[[] for i in xrange(num_popln)]
num_spks_per_cell_per_interbin=[[] for i in xrange(num_popln)]
num_spks_per_cell_per_interbin_no_trans=[[] for i in xrange(num_popln)]
avg_num_spks_per_cell_inter=[[] for i in xrange(num_popln)]
std_num_spks_per_cell_inter=[[] for i in xrange(num_popln)]
tot_num_spks_per_interbin=[[] for i in xrange(num_popln)]
tot_num_spks_per_interbin_no_trans=[[] for i in xrange(num_popln)]
avg_tot_num_spks_inter=[[] for i in xrange(num_popln)]
std_tot_num_spks_inter=[[] for i in xrange(num_popln)]
avg_interburst_freq_per_bin_tot=[[] for i in xrange(num_popln)]
std_interburst_freq_per_bin_tot=[[] for i in xrange(num_popln)]
avg_interburst_freq_per_cell_tot=[[] for i in xrange(num_popln)]
std_interburst_freq_per_cell_tot=[[] for i in xrange(num_popln)]
avg_interburst_freq_per_bin=[[] for i in xrange(num_popln)]
std_interburst_freq_per_bin=[[] for i in xrange(num_popln)]
avg_interBF_per_bin_list=[[] for i in xrange(num_popln)]
std_interBF_per_bin_list=[[] for i in xrange(num_popln)]
avg_interBF_per_cell_list=[[] for i in xrange(num_popln)]
std_interBF_per_cell_list=[[] for i in xrange(num_popln)]

###########################################################
###### determine stats for PYR cells (the number of spikes per burst, number of cells spiking, etc).
if (numbins-bin_ind_no_trans-1)>1:  #must have at least two bins following transient to do stats (since last one will not count and need at least one bin)
    for i in xrange(num_popln):
        [indices[i],num_cells_w_FBF[i],fullburst_freq_nonzero[i],fullburst_freq_nonzero_no_trans[i],
                num_cells_spk_per_bin[i],num_cells_spk_per_bin_no_trans[i],avg_num_cells_spk[i],
                std_num_cells_spk[i],num_spks_per_cell_per_bin[i],num_spks_per_cell_per_bin_no_trans[i],
                avg_num_spks_per_cell[i],std_num_spks_per_cell[i],tot_num_spks_per_bin[i],
                tot_num_spks_per_bin_no_trans[i],avg_tot_num_spks[i],std_tot_num_spks[i],avg_fullburst_freq_per_bin_tot[i],
                std_fullburst_freq_per_bin_tot[i],avg_fullburst_freq_per_cell_tot[i],std_fullburst_freq_per_cell_tot[i],
                avg_fullburst_freq_per_bin[i],std_fullburst_freq_per_bin[i],
                avg_FBF_per_bin_list[i],std_FBF_per_bin_list[i],avg_FBF_per_cell_list[i],std_FBF_per_cell_list[i],
                intraindices[i],num_cells_w_intraBF[i],intraburst_freq_nonzero[i],intraburst_freq_nonzero_no_trans[i], 
                num_cells_spk_per_intrabin[i],num_cells_spk_per_intrabin_no_trans[i],avg_num_cells_spk_intra[i],
                std_num_cells_spk_intra[i],num_spks_per_cell_per_intrabin[i],num_spks_per_cell_per_intrabin_no_trans[i],
                avg_num_spks_per_cell_intra[i],std_num_spks_per_cell_intra[i],tot_num_spks_per_intrabin[i],
                tot_num_spks_per_intrabin_no_trans[i],avg_tot_num_spks_intra[i],std_tot_num_spks_intra[i],avg_intraburst_freq_per_bin_tot[i],
                std_intraburst_freq_per_bin_tot[i],avg_intraburst_freq_per_cell_tot[i],std_intraburst_freq_per_cell_tot[i],
                avg_intraburst_freq_per_bin[i],std_intraburst_freq_per_bin[i],
                avg_intraBF_per_bin_list[i],std_intraBF_per_bin_list[i],avg_intraBF_per_cell_list[i],std_intraBF_per_cell_list[i],
                interindices[i],num_cells_w_interBF[i],interburst_freq_nonzero[i],interburst_freq_nonzero_no_trans[i], 
                num_cells_spk_per_interbin[i],num_cells_spk_per_interbin_no_trans[i],avg_num_cells_spk_inter[i],
                std_num_cells_spk_inter[i],num_spks_per_cell_per_interbin[i],num_spks_per_cell_per_interbin_no_trans[i],
                avg_num_spks_per_cell_inter[i],std_num_spks_per_cell_inter[i],tot_num_spks_per_interbin[i],
                tot_num_spks_per_interbin_no_trans[i],avg_tot_num_spks_inter[i],std_tot_num_spks_inter[i],avg_interburst_freq_per_bin_tot[i],
                std_interburst_freq_per_bin_tot[i],avg_interburst_freq_per_cell_tot[i],std_interburst_freq_per_cell_tot[i],
                avg_interburst_freq_per_bin[i],std_interburst_freq_per_bin[i],
                avg_interBF_per_bin_list[i],std_interBF_per_bin_list[i],avg_interBF_per_cell_list[i],
                std_interBF_per_cell_list[i]]= bursting_stats(bin_s,bin_ms,binwidth_ms,numbins,intraburst_bins/ms,interburst_bins/ms,
                                                              num_intraburst_bins,num_interburst_bins,
                                                              bin_ind_no_trans,intrabin_ind_no_trans,interbin_ind_no_trans,
                                                              N_list[i],M_t_list[i]/ms,M_i_list[i],transient)

##    [indices,num_cells_w_FBF,fullburst_freq_nonzero,fullburst_freq_nonzero_no_trans,
##        num_cells_spk_per_bin,num_cells_spk_per_bin_no_trans,avg_num_cells_spk,
##        std_num_cells_spk,num_spks_per_cell_per_bin,num_spks_per_cell_per_bin_no_trans,
##        avg_num_spks_per_cell,std_num_spks_per_cell,tot_num_spks_per_bin,
##        tot_num_spks_per_bin_no_trans,avg_tot_num_spks,std_tot_num_spks,avg_fullburst_freq_per_bin_tot,
##        std_fullburst_freq_per_bin_tot,avg_fullburst_freq_per_cell_tot,std_fullburst_freq_per_cell_tot,
##        avg_fullburst_freq_per_bin,std_fullburst_freq_per_bin,
##        avg_FBF_per_bin_list,std_FBF_per_bin_list,avg_FBF_per_cell_list,std_FBF_per_cell_list,
##        intraindices,num_cells_w_intraBF,intraburst_freq_nonzero,intraburst_freq_nonzero_no_trans, 
##        num_cells_spk_per_intrabin,num_cells_spk_per_intrabin_no_trans,avg_num_cells_spk_intra,
##        std_num_cells_spk_intra,num_spks_per_cell_per_intrabin,num_spks_per_cell_per_intrabin_no_trans,
##        avg_num_spks_per_cell_intra,std_num_spks_per_cell_intra,tot_num_spks_per_intrabin,
##        tot_num_spks_per_intrabin_no_trans,avg_tot_num_spks_intra,std_tot_num_spks_intra,avg_intraburst_freq_per_bin_tot,
##        std_intraburst_freq_per_bin_tot,avg_intraburst_freq_per_cell_tot,std_intraburst_freq_per_cell_tot,
##        avg_intraburst_freq_per_bin,std_intraburst_freq_per_bin,
##        avg_intraBF_per_bin_list,std_intraBF_per_bin_list,avg_intraBF_per_cell_list,std_intraBF_per_cell_list,
##        interindices,num_cells_w_interBF,interburst_freq_nonzero,interburst_freq_nonzero_no_trans, 
##        num_cells_spk_per_interbin,num_cells_spk_per_interbin_no_trans,avg_num_cells_spk_inter,
##        std_num_cells_spk_inter,num_spks_per_cell_per_interbin,num_spks_per_cell_per_interbin_no_trans,
##        avg_num_spks_per_cell_inter,std_num_spks_per_cell_inter,tot_num_spks_per_interbin,
##        tot_num_spks_per_interbin_no_trans,avg_tot_num_spks_inter,std_tot_num_spks_inter,avg_interburst_freq_per_bin_tot,
##        std_interburst_freq_per_bin_tot,avg_interburst_freq_per_cell_tot,std_interburst_freq_per_cell_tot,
##        avg_interburst_freq_per_bin,std_interburst_freq_per_bin,
##        avg_interBF_per_bin_list,std_interBF_per_bin_list,avg_interBF_per_cell_list,
##        std_interBF_per_cell_list] = bursting_stats(bin_s,bin_ms,binwidth_ms,numbins,intraburst_bins/ms,interburst_bins/ms,
##                                                      num_intraburst_bins,num_interburst_bins,
##                                                      bin_ind_no_trans,intrabin_ind_no_trans,interbin_ind_no_trans,
##                                                      N_pyr,Spiketimes_PYR.t/ms,Spiketimes_PYR.i,transient)

    print("Writing stats")
    write_stats(bin_ms,binpts,num_popln,avg_num_cells_spk,std_num_cells_spk,
                    avg_num_spks_per_cell,std_num_spks_per_cell,avg_tot_num_spks,std_tot_num_spks,
                    avg_fullburst_freq_per_bin_tot,std_fullburst_freq_per_bin_tot,
                    totspkhist_list,totspkdist_smooth,num_cells_spk_per_bin,num_spks_per_cell_per_bin,
                    tot_num_spks_per_bin,avg_FBF_per_bin_list,std_FBF_per_bin_list,num_cells_w_FBF,
                    avg_fullburst_freq_per_cell_tot,std_fullburst_freq_per_cell_tot,
                    indices,avg_FBF_per_cell_list,std_FBF_per_cell_list,bin_ind_no_trans,
                    numbins,fullburst_freq_nonzero,
                    intraburst_bins_ms,avg_num_cells_spk_intra,std_num_cells_spk_intra,
                    avg_num_spks_per_cell_intra,std_num_spks_per_cell_intra,avg_tot_num_spks_intra,std_tot_num_spks_intra,
                    avg_intraburst_freq_per_bin_tot,std_intraburst_freq_per_bin_tot,
                    num_cells_spk_per_intrabin,num_spks_per_cell_per_intrabin,
                    tot_num_spks_per_intrabin,avg_intraBF_per_bin_list,std_intraBF_per_bin_list,num_cells_w_intraBF,
                    avg_intraburst_freq_per_cell_tot,std_intraburst_freq_per_cell_tot,
                    intraindices,avg_intraBF_per_cell_list,std_intraBF_per_cell_list,intrabin_ind_no_trans,
                    num_intraburst_bins,intraburst_freq_nonzero,
                    interburst_bins_ms,avg_num_cells_spk_inter,std_num_cells_spk_inter,
                    avg_num_spks_per_cell_inter,std_num_spks_per_cell_inter,avg_tot_num_spks_inter,std_tot_num_spks_inter,
                    avg_interburst_freq_per_bin_tot,std_interburst_freq_per_bin_tot,
                    num_cells_spk_per_interbin,num_spks_per_cell_per_interbin,
                    tot_num_spks_per_interbin,avg_interBF_per_bin_list,std_interBF_per_bin_list,num_cells_w_interBF,
                    avg_interburst_freq_per_cell_tot,std_interburst_freq_per_cell_tot,
                    interindices,avg_interBF_per_cell_list,std_interBF_per_cell_list,interbin_ind_no_trans,
                    num_interburst_bins,interburst_freq_nonzero,filename,foldername)
 
    #Plots
    make_stats_plots(num_popln,N_list,ctr_of_bin_ms,ctr_of_intrabin_ms,ctr_of_interbin_ms,transient,
                         num_cells_spk_per_bin,avg_num_cells_spk,std_num_cells_spk,
                         num_spks_per_cell_per_bin,avg_num_spks_per_cell,std_num_spks_per_cell,
                         tot_num_spks_per_bin,avg_tot_num_spks,std_tot_num_spks,
                         avg_fullburst_freq_per_bin,std_fullburst_freq_per_bin,
                         avg_FBF_per_bin_list,
                         avg_fullburst_freq_per_bin_tot,std_fullburst_freq_per_bin_tot,
                         avg_FBF_per_cell_list,std_FBF_per_cell_list,
                         avg_fullburst_freq_per_cell_tot,std_fullburst_freq_per_cell_tot,
                         num_cells_spk_per_intrabin,avg_num_cells_spk_intra,std_num_cells_spk_intra,
                         num_spks_per_cell_per_intrabin,avg_num_spks_per_cell_intra,std_num_spks_per_cell_intra,
                         tot_num_spks_per_intrabin,avg_tot_num_spks_intra,std_tot_num_spks_intra,
                         avg_intraburst_freq_per_bin,std_intraburst_freq_per_bin,
                         avg_intraBF_per_bin_list,
                         avg_intraburst_freq_per_bin_tot,std_intraburst_freq_per_bin_tot,
                         avg_intraBF_per_cell_list,std_intraBF_per_cell_list,
                         avg_intraburst_freq_per_cell_tot,std_intraburst_freq_per_cell_tot,
                         num_cells_spk_per_interbin,avg_num_cells_spk_inter,std_num_cells_spk_inter,
                         num_spks_per_cell_per_interbin,avg_num_spks_per_cell_inter,std_num_spks_per_cell_inter,
                         tot_num_spks_per_interbin,avg_tot_num_spks_inter,std_tot_num_spks_inter,
                         avg_interburst_freq_per_bin,std_interburst_freq_per_bin,
                         avg_interBF_per_bin_list,
                         avg_interburst_freq_per_bin_tot,std_interburst_freq_per_bin_tot,
                         avg_interBF_per_cell_list,std_interBF_per_cell_list,
                         avg_interburst_freq_per_cell_tot,std_interburst_freq_per_cell_tot,
                         celltype,filename,foldername)


#Save power spectrum information
fps=open(foldername+"/"+filename+"_power_spec.txt","w")
fps.write(str(max_freq) + " " + str(max_pwr))
fps.close()

plot(freq[2:length/2],network_fft[2:length/2])
xlim([0,80])
xlabel("Frequency (Hz)")
ylabel("Power")
savefig(foldername+"/"+filename+"_power_spec.png")
close()

#Num of models
fMods = open(foldername+"/"+filename+"_models.txt","w")
fMods.write(str(lenOfAPS))
fMods.close()

#Save rasterplots
plot(Spiketimes_PYR.t/ms,Spiketimes_PYR.i,'.')    
xlabel("Time (ms)")
ylabel("Cell Number")
savefig(foldername+"/"+filename+"_PYR_raster.png")
close()

plot(Spiketimes_PV.t/ms,Spiketimes_PV.i,'.')
xlabel("Time (ms)")
ylabel("Cell Number")
savefig(foldername+'/'+filename+'_PV_raster.png')
close()

plot(Spiketimes_PYR.t[Spiketimes_PYR.t >= duration - (1 * second)]/ms,Spiketimes_PYR.i[Spiketimes_PYR.t >= duration - (1 * second)],'.')    
xlabel("Time (ms)")
ylabel("Cell Number")
savefig(foldername+"/"+filename+"_PYR_raster_end.png")
close()

plot(Spiketimes_PV.t[Spiketimes_PV.t >= duration - (1 * second)]/ms,Spiketimes_PV.i[Spiketimes_PV.t >= duration - (1 * second)],'.')    
xlabel("Time (ms)")
ylabel("Cell Number")
savefig(foldername+"/"+filename+"_PV_raster_end.png")
close()


#Peakfinding and stats
all_PV_E_peaks = list()
all_PV_I_peaks = list()

all_PYR_E_peaks = list()
all_PYR_I_peaks = list()

for i in range(N_pv_mon):
    PV_E = abs((PV_Esyn.Isyn_pyr_pv[i])/pA)[PV_Esyn.t > transient * ms]
    PV_I = abs((PV_Isyn.Isyn_pv_pv[i])/pA)[PV_Esyn.t > transient * ms]
    peaks_PV_E = find_peaks_cwt(PV_E,arange(30,150))
    peaks_PV_I = find_peaks_cwt(PV_I,arange(30,150))
    all_PV_E_peaks.extend(PV_E[peaks_PV_E])
    all_PV_I_peaks.extend(PV_I[peaks_PV_I])

for i in range(N_pyr_mon):
    PYR_E = abs((PYR_Esyn.Isyn_pyr_pyr[i])/pA)[PV_Esyn.t > transient * ms]
    PYR_I = abs((PYR_Isyn.Isyn_pv_pyr[i])/pA)[PV_Esyn.t > transient * ms]
    peaks_PYR_E = find_peaks_cwt(PYR_E,arange(30,150))
    peaks_PYR_I = find_peaks_cwt(PYR_I,arange(30,150))
    all_PYR_E_peaks.extend(PYR_E[peaks_PYR_E])
    all_PYR_I_peaks.extend(PYR_I[peaks_PYR_I])
    


#Stats on PSC
avg_EPSC_PV = mean(all_PV_E_peaks)
std_EPSC_PV = std(all_PV_E_peaks)

avg_IPSC_PV = mean(all_PV_I_peaks)
std_IPSC_PV = std(all_PV_I_peaks)

avg_EPSC_PYR = mean(all_PYR_E_peaks)
std_EPSC_PYR = std(all_PYR_E_peaks)

avg_IPSC_PYR = mean(all_PYR_I_peaks)
std_IPSC_PYR = std(all_PYR_I_peaks)

fPSC_averages=open(foldername+"/"+filename+"_PSC_averages.txt","w")

#E/I ratio for PYR
if (logical_and(avg_IPSC_PYR != 0.0,avg_IPSC_PYR == avg_IPSC_PYR)):
    fPSC_averages.write(str(avg_EPSC_PYR/avg_IPSC_PYR))
else:
    fPSC_averages.write(str('Inv'))
fPSC_averages.write(str('\n'))
#E/I ratio of PV
if (logical_and(avg_IPSC_PV != 0.0,avg_IPSC_PV == avg_IPSC_PV)):
    fPSC_averages.write(str(avg_EPSC_PV/avg_IPSC_PV))
else:
    fPSC_averages.write(str('Inv'))
fPSC_averages.write(str('\n'))

#Averages of PYR and PV for EPSC and IPSC
fPSC_averages.write(str(avg_EPSC_PYR))
fPSC_averages.write('+/-')
fPSC_averages.write(str(std_EPSC_PYR))
fPSC_averages.write("\n")

fPSC_averages.write(str(avg_IPSC_PYR))
fPSC_averages.write('+/-')
fPSC_averages.write(str(std_IPSC_PYR))
fPSC_averages.write("\n")

fPSC_averages.write(str(avg_EPSC_PV))
fPSC_averages.write('+/-')
fPSC_averages.write(str(std_EPSC_PV))
fPSC_averages.write("\n")

fPSC_averages.write(str(avg_IPSC_PV))
fPSC_averages.write('+/-')
fPSC_averages.write(str(std_IPSC_PV))
fPSC_averages.write("\n")

fPSC_averages.close()

#save('./Output/Isyn_pyr_pyr',PYR_Esyn.Isyn_pyr_pyr/pA)

print("Done!")


