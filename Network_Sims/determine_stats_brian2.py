
'''
Created on 2014-05-02

@author: kferguso
Edited by Anton Lunyov on 1 Aug, 2017 to work with Brian2
'''
from brian2 import *

from scipy.signal import filtfilt,butter
from numpy.ma.core import zeros,ceil,logical_and, where, reshape, floor, argmax, mean, std, exp, nonzero
import bisect
from pylab import average
from matplotlib.pyplot import figure, plot, ylim, xlim, ylabel, xlabel, title, show, savefig, close, suptitle, errorbar, subplot,tight_layout
from matplotlib.ticker import MultipleLocator

##### determine stats: the number of spikes per burst, number of cells spiking, etc.
##collect all stats, but find averages, etc by eliminating transient
def bursting_stats(bin_s,bin_ms,binwidth_ms,numbins,intraburst_bins,interburst_bins,num_intraburst_bins,num_interburst_bins,
				   bin_ind_no_trans,intrabin_ind_no_trans,interbin_ind_no_trans,N,M_t,M_i,transient):
	numspks=zeros((numbins,N))
	fullburst_freq=zeros((numbins,N))
	
	#since distn always starts on a decrease, you either have the same number of interbursts as intra (if end on a dec), or one more inter than intra (if end on an inc) 
	numspks_intra=zeros((num_intraburst_bins,N))
	numspks_inter=zeros((num_interburst_bins,N))
	intraburst_freq=zeros((num_intraburst_bins,N))
	interburst_freq=zeros((num_interburst_bins,N))
	   
	for i in xrange(numbins):
		step_start=bin_s[i] * 1000
		step_end=bin_s[i+1] * 1000    
		
		for j in xrange(N):
			spks=where(logical_and(logical_and(M_t>step_start, M_t<step_end),M_i == j))
			numspks[i,j] = sum(logical_and(logical_and(M_t>step_start, M_t<step_end),M_i == j))
			
			#find fullburst freq
			if numspks[i,j]>1:
				spktms=M_t[spks]
				fullburst_freq[i,j]=1/((spktms[1:]-spktms[:-1]).mean(axis=0))
	
	for i in xrange(num_intraburst_bins): 
		step_start_intra=intraburst_bins[2*i] * second
		step_start_inter=interburst_bins[2*i] * second
		step_end_intra=intraburst_bins[2*i+1] * second
		step_end_inter=interburst_bins[2*i+1] * second
		
		for j in xrange(N):
			spks_intra=where(logical_and(logical_and(M_t>step_start_intra, M_t<step_end_intra),M_i == j))
			spks_inter=where(logical_and(logical_and(M_t>step_start_inter, M_t<step_end_inter),M_i == j))
			numspks_intra[i,j]=sum(logical_and(logical_and(M_t>step_start_intra, M_t<step_end_intra),M_i == j))
			numspks_inter[i,j]=sum(logical_and(logical_and(M_t>step_start_inter, M_t<step_end_inter),M_i == j))
			
			#find intraburst and interburst freq
			if numspks_intra[i,j]>1:
				spktms_intra=M_t[spks_intra]
				intraburst_freq[i,j]=1/((spktms_intra[1:]-spktms_intra[:-1]).mean(axis=0))
			if numspks_inter[i,j]>1:
				spktms_inter=M_t[spks_inter]
				interburst_freq[i,j]=1/((spktms_inter[1:]-spktms_inter[:-1]).mean(axis=0))    
	
	#just for last value (will have one extra if more inter than intra)
	if (num_interburst_bins>num_intraburst_bins):
		step_start_inter=interburst_bins[-2] * second
		step_end_inter=interburst_bins[-1] * second
		
		for j in xrange(N):
			spks_inter=where(logical_and(logical_and(M_t>step_start_inter, M_t<step_end_inter),M_i == j))
			numspks_inter[-1,j]=sum(logical_and(logical_and(M_t>step_start_inter, M_t<step_end_inter),M_i == j))
			
			if numspks_inter[-1,j]>1:
				spktms_inter=M_t[spks_inter]
				interburst_freq[-1,j]=1/((spktms_inter[1:]-spktms_inter[:-1]).mean(axis=0))
		
		   
	num_cells_spk_per_bin=(numspks!=0).sum(1)
	num_cells_spk_per_intrabin=(numspks_intra!=0).sum(1)
	num_cells_spk_per_interbin=(numspks_inter!=0).sum(1)
	
	#print(num_cells_spk_per_bin)

	
#    #find first index after transient
#    bin_ind_no_trans=bisect.bisect(bin_ms, transient)
#    intrabin_ind_no_trans=bisect.bisect(intraburst_bins, transient/1000)  #transient to seconds      
#    if intrabin_ind_no_trans % 2 != 0:   #index must be even since format is ind0=start_bin, ind1=end_bin, ind2=start_bin, .... .
#        intrabin_ind_no_trans += 1
#    interbin_ind_no_trans=bisect.bisect(interburst_bins, transient/1000)
#    if interbin_ind_no_trans % 2 != 0:
#        interbin_ind_no_trans += 1
	
	#find averages and stds after transient
	#for full burst
	num_cells_spk_per_bin_no_trans=num_cells_spk_per_bin[bin_ind_no_trans:-2]  #ignore transient and last bin
	avg_num_cells_spk=num_cells_spk_per_bin_no_trans.mean(axis=0)
	std_num_cells_spk=num_cells_spk_per_bin_no_trans.std(axis=0)
	
	#for intra burst
	num_cells_spk_per_intrabin_no_trans=num_cells_spk_per_intrabin[intrabin_ind_no_trans:-1]  #ignore transient and last bin
	avg_num_cells_spk_intra=num_cells_spk_per_intrabin_no_trans.mean(axis=0)
	std_num_cells_spk_intra=num_cells_spk_per_intrabin_no_trans.std(axis=0)
	
	#for inter burst
	num_cells_spk_per_interbin_no_trans=num_cells_spk_per_interbin[interbin_ind_no_trans:-1]  #ignore transient and last bin
	avg_num_cells_spk_inter=num_cells_spk_per_interbin_no_trans.mean(axis=0)
	std_num_cells_spk_inter=num_cells_spk_per_interbin_no_trans.std(axis=0)
	
	#for full burst
	num_spks_per_cell_per_bin=numspks.mean(axis=1)
	num_spks_per_cell_per_bin_no_trans=num_spks_per_cell_per_bin[bin_ind_no_trans:-2]
	avg_num_spks_per_cell=num_spks_per_cell_per_bin_no_trans.mean(axis=0)
	std_num_spks_per_cell=num_spks_per_cell_per_bin_no_trans.std(axis=0)
	
	#for intra burst
	num_spks_per_cell_per_intrabin=numspks_intra.mean(axis=1)
	num_spks_per_cell_per_intrabin_no_trans=num_spks_per_cell_per_intrabin[intrabin_ind_no_trans:-1]
	avg_num_spks_per_cell_intra=num_spks_per_cell_per_intrabin_no_trans.mean(axis=0)
	std_num_spks_per_cell_intra=num_spks_per_cell_per_intrabin_no_trans.std(axis=0)
	
	#for inter burst
	num_spks_per_cell_per_interbin=numspks_inter.mean(axis=1)
	num_spks_per_cell_per_interbin_no_trans=num_spks_per_cell_per_interbin[interbin_ind_no_trans:-1]
	avg_num_spks_per_cell_inter=num_spks_per_cell_per_interbin_no_trans.mean(axis=0)
	std_num_spks_per_cell_inter=num_spks_per_cell_per_interbin_no_trans.std(axis=0)
	
	#for full burst
	tot_num_spks_per_bin=numspks.sum(1)
	tot_num_spks_per_bin_no_trans=tot_num_spks_per_bin[bin_ind_no_trans:-2]
	avg_tot_num_spks=tot_num_spks_per_bin_no_trans.mean(axis=0)
	std_tot_num_spks=tot_num_spks_per_bin_no_trans.std(axis=0)
	
	#for intra burst
	tot_num_spks_per_intrabin=numspks_intra.sum(1)
	tot_num_spks_per_intrabin_no_trans=tot_num_spks_per_intrabin[intrabin_ind_no_trans:-1]
	avg_tot_num_spks_intra=tot_num_spks_per_intrabin_no_trans.mean(axis=0)
	std_tot_num_spks_intra=tot_num_spks_per_intrabin_no_trans.std(axis=0)
	
	#for inter burst
	tot_num_spks_per_interbin=numspks_inter.sum(1)
	tot_num_spks_per_interbin_no_trans=tot_num_spks_per_interbin[interbin_ind_no_trans:-1]
	avg_tot_num_spks_inter=tot_num_spks_per_interbin_no_trans.mean(axis=0)
	std_tot_num_spks_inter=tot_num_spks_per_interbin_no_trans.std(axis=0)
	
	#for full burst
	avg_fullburst_freq_per_bin=zeros((numbins,1))   #can do these stats on all bins and ignore trans for total average 
	std_fullburst_freq_per_bin=zeros((numbins,1))    
	avg_fullburst_freq_per_cell=zeros((N,1))      #must ignore trans immediately here... 
	std_fullburst_freq_per_cell=zeros((N,1))
	#for intra burst
	avg_intraburst_freq_per_bin=zeros((num_intraburst_bins,1))   #can do these stats on all bins and ignore trans for total average 
	std_intraburst_freq_per_bin=zeros((num_intraburst_bins,1))    
	avg_intraburst_freq_per_cell=zeros((N,1))      #must ignore trans immediately here... 
	std_intraburst_freq_per_cell=zeros((N,1))
	# for inter burst
	avg_interburst_freq_per_bin=zeros((num_interburst_bins,1))   #can do these stats on all bins and ignore trans for total average 
	std_interburst_freq_per_bin=zeros((num_interburst_bins,1))    
	avg_interburst_freq_per_cell=zeros((N,1))      #must ignore trans immediately here... 
	std_interburst_freq_per_cell=zeros((N,1))
	
	#for full burst
	FBF_ind=nonzero(fullburst_freq)
	fullburst_freq_no_trans=fullburst_freq[bin_ind_no_trans:-2][:]
	FBF_ind_no_trans=nonzero(fullburst_freq_no_trans)
	uniq_ind_per_bin=sorted(set(FBF_ind[0]))    #rows = bins
	uniq_ind_per_cell=sorted(set(FBF_ind_no_trans[1]))    #colms = cell
	#for intraburst
	intraBF_ind=nonzero(intraburst_freq)
	intraburst_freq_no_trans=intraburst_freq[intrabin_ind_no_trans:-1][:]
	intraBF_ind_no_trans=nonzero(intraburst_freq_no_trans)
	uniq_ind_per_intrabin=sorted(set(intraBF_ind[0]))    #rows = bins
	uniq_ind_per_cell_intra=sorted(set(intraBF_ind_no_trans[1]))    #colms = cell
	# for interburst
	interBF_ind=nonzero(interburst_freq)
	interburst_freq_no_trans=interburst_freq[interbin_ind_no_trans:-1][:]
	interBF_ind_no_trans=nonzero(interburst_freq_no_trans)
	uniq_ind_per_interbin=sorted(set(interBF_ind[0]))    #rows = bins
	uniq_ind_per_cell_inter=sorted(set(interBF_ind_no_trans[1]))    #colms = cell
	
	for i in uniq_ind_per_bin:
		nonzero_row_ind=[j for j, x in enumerate(FBF_ind[0]) if x==i]
		FBF_nonzero_row=[]
		for k in xrange(len(nonzero_row_ind)):
			FBF_nonzero_row.append(fullburst_freq[FBF_ind[0][nonzero_row_ind[k]]][FBF_ind[1][nonzero_row_ind[k]]])
		avg_fullburst_freq_per_bin[i]=average(FBF_nonzero_row)
		std_fullburst_freq_per_bin[i]=std(FBF_nonzero_row)
	for i in uniq_ind_per_intrabin:
		nonzero_row_intraind=[j for j, x in enumerate(intraBF_ind[0]) if x==i]
		intraBF_nonzero_row=[]
		for k in xrange(len(nonzero_row_intraind)):
			intraBF_nonzero_row.append(intraburst_freq[intraBF_ind[0][nonzero_row_intraind[k]]][intraBF_ind[1][nonzero_row_intraind[k]]])
		avg_intraburst_freq_per_bin[i]=average(intraBF_nonzero_row)
		std_intraburst_freq_per_bin[i]=std(intraBF_nonzero_row)
	for i in uniq_ind_per_interbin:
		nonzero_row_interind=[j for j, x in enumerate(interBF_ind[0]) if x==i]
		interBF_nonzero_row=[]
		for k in xrange(len(nonzero_row_interind)):
			interBF_nonzero_row.append(interburst_freq[interBF_ind[0][nonzero_row_interind[k]]][interBF_ind[1][nonzero_row_interind[k]]])
		avg_interburst_freq_per_bin[i]=average(interBF_nonzero_row)
		std_interburst_freq_per_bin[i]=std(interBF_nonzero_row)

	
	for i in uniq_ind_per_cell:
		nonzero_col_ind=[j for j, x in enumerate(FBF_ind_no_trans[1]) if x==i]
		FBF_nonzero_col_no_trans=[]
		for k in xrange(len(nonzero_col_ind)):
			FBF_nonzero_col_no_trans.append(fullburst_freq_no_trans[FBF_ind_no_trans[0][nonzero_col_ind[k]]][FBF_ind_no_trans[1][nonzero_col_ind[k]]])
		avg_fullburst_freq_per_cell[i]=average(FBF_nonzero_col_no_trans)
		std_fullburst_freq_per_cell[i]=std(FBF_nonzero_col_no_trans) 
	for i in uniq_ind_per_cell_intra:
		nonzero_col_intraind=[j for j, x in enumerate(intraBF_ind_no_trans[1]) if x==i]
		intraBF_nonzero_col_no_trans=[]
		for k in xrange(len(nonzero_col_intraind)):
			intraBF_nonzero_col_no_trans.append(intraburst_freq_no_trans[intraBF_ind_no_trans[0][nonzero_col_intraind[k]]][intraBF_ind_no_trans[1][nonzero_col_intraind[k]]])
		avg_intraburst_freq_per_cell[i]=average(intraBF_nonzero_col_no_trans)
		std_intraburst_freq_per_cell[i]=std(intraBF_nonzero_col_no_trans)
	for i in uniq_ind_per_cell_inter:
		nonzero_col_interind=[j for j, x in enumerate(interBF_ind_no_trans[1]) if x==i]
		interBF_nonzero_col_no_trans=[]
		for k in xrange(len(nonzero_col_interind)):
			interBF_nonzero_col_no_trans.append(interburst_freq_no_trans[interBF_ind_no_trans[0][nonzero_col_interind[k]]][interBF_ind_no_trans[1][nonzero_col_interind[k]]])
		avg_interburst_freq_per_cell[i]=average(interBF_nonzero_col_no_trans)
		std_interburst_freq_per_cell[i]=std(interBF_nonzero_col_no_trans)
	
	#Find avg FBF. Only find average of bins that have an FBF (i.e. eliminate non-zero bins).  Here we have found the average of the average FBF/bin.  That is, we have represented each bin equally, rather than each set of spikes....
	if len(avg_fullburst_freq_per_bin)>(bin_ind_no_trans+2):  #can we eliminate transient and last bin? 
		avg_fullburst_freq_per_bin_no_trans=avg_fullburst_freq_per_bin[bin_ind_no_trans:-2]   #must eliminate trans here since didn't before
		if max(avg_fullburst_freq_per_bin_no_trans)>0:    
			avg_fullburst_freq_per_bin_tot=(avg_fullburst_freq_per_bin_no_trans[nonzero(avg_fullburst_freq_per_bin_no_trans)]).mean(axis=0)
			std_fullburst_freq_per_bin_tot=(avg_fullburst_freq_per_bin_no_trans[nonzero(avg_fullburst_freq_per_bin_no_trans)]).std(axis=0)
		else:
			avg_fullburst_freq_per_bin_tot=0
			std_fullburst_freq_per_bin_tot=0
	else:
		avg_fullburst_freq_per_bin_no_trans=0
		avg_fullburst_freq_per_bin_tot=0
		std_fullburst_freq_per_bin_tot=0
			
	if max(avg_fullburst_freq_per_cell)>0:        
		avg_fullburst_freq_per_cell_tot=(avg_fullburst_freq_per_cell[nonzero(avg_fullburst_freq_per_cell)]).mean(axis=0)
		std_fullburst_freq_per_cell_tot=(avg_fullburst_freq_per_cell[nonzero(avg_fullburst_freq_per_cell)]).std(axis=0)
	else:
		avg_fullburst_freq_per_cell_tot=0
		std_fullburst_freq_per_cell_tot=0
		
	#intra
	if len(avg_intraburst_freq_per_bin)>(intrabin_ind_no_trans+2):  #can we eliminate transient and last bin? 
		avg_intraburst_freq_per_bin_no_trans=avg_intraburst_freq_per_bin[intrabin_ind_no_trans:-1]   #must eliminate trans here since didn't before
		if max(avg_intraburst_freq_per_bin_no_trans)>0:    
			avg_intraburst_freq_per_bin_tot=(avg_intraburst_freq_per_bin_no_trans[nonzero(avg_intraburst_freq_per_bin_no_trans)]).mean(axis=0)
			std_intraburst_freq_per_bin_tot=(avg_intraburst_freq_per_bin_no_trans[nonzero(avg_intraburst_freq_per_bin_no_trans)]).std(axis=0)
		else:
			avg_intraburst_freq_per_bin_tot=0
			std_intraburst_freq_per_bin_tot=0
	else:
		avg_intraburst_freq_per_bin_no_trans=0
		avg_intraburst_freq_per_bin_tot=0
		std_intraburst_freq_per_bin_tot=0
		
	if max(avg_intraburst_freq_per_cell)>0:    
		avg_intraburst_freq_per_cell_tot=(avg_intraburst_freq_per_cell[nonzero(avg_intraburst_freq_per_cell)]).mean(axis=0)
		std_intraburst_freq_per_cell_tot=(avg_intraburst_freq_per_cell[nonzero(avg_intraburst_freq_per_cell)]).std(axis=0)
	else:
		avg_intraburst_freq_per_cell_tot=0
		std_intraburst_freq_per_cell_tot=0

  
	#inter
	if len(avg_interburst_freq_per_bin)>(interbin_ind_no_trans+2):  #can we eliminate transient and last bin? 
		avg_interburst_freq_per_bin_no_trans=avg_interburst_freq_per_bin[interbin_ind_no_trans:-1]   #must eliminate trans here since didn't before
		if max(avg_interburst_freq_per_bin_no_trans)>0:    
			avg_interburst_freq_per_bin_tot=(avg_interburst_freq_per_bin_no_trans[nonzero(avg_interburst_freq_per_bin_no_trans)]).mean(axis=0)
			std_interburst_freq_per_bin_tot=(avg_interburst_freq_per_bin_no_trans[nonzero(avg_interburst_freq_per_bin_no_trans)]).std(axis=0)
		else:
			avg_interburst_freq_per_bin_tot=0
			std_interburst_freq_per_bin_tot=0
	else:
		avg_interburst_freq_per_bin_no_trans=0
		avg_interburst_freq_per_bin_tot=0
		std_interburst_freq_per_bin_tot=0
	
	if max(avg_interburst_freq_per_cell)>0:    
		avg_interburst_freq_per_cell_tot=(avg_interburst_freq_per_cell[nonzero(avg_interburst_freq_per_cell)]).mean(axis=0)
		std_interburst_freq_per_cell_tot=(avg_interburst_freq_per_cell[nonzero(avg_interburst_freq_per_cell)]).std(axis=0)
	else:
		avg_interburst_freq_per_cell_tot=0
		std_interburst_freq_per_cell_tot=0
	
	#convert FBF to list for writing and plotting
	avg_FBF_per_bin_list=[val for subl in avg_fullburst_freq_per_bin for val in subl]
	std_FBF_per_bin_list=[val for subl in std_fullburst_freq_per_bin for val in subl]
	avg_FBF_per_cell_list=[val for subl in avg_fullburst_freq_per_cell for val in subl]
	std_FBF_per_cell_list=[val for subl in std_fullburst_freq_per_cell for val in subl]
  
	indices_full=nonzero(fullburst_freq)
	indices=sorted(set(indices_full[1]))
	fullburst_freq_nonzero=fullburst_freq[:,indices] #[1]]
	fullburst_freq_nonzero_no_trans=fullburst_freq_nonzero[bin_ind_no_trans:-2][:]
	num_cells_w_FBF=len(indices) #[1]) #len(fullburst_freq_nonzero_no_trans)
	
	avg_intraBF_per_bin_list=[val for subl in avg_intraburst_freq_per_bin for val in subl]
	std_intraBF_per_bin_list=[val for subl in std_intraburst_freq_per_bin for val in subl]
	avg_intraBF_per_cell_list=[val for subl in avg_intraburst_freq_per_cell for val in subl]
	std_intraBF_per_cell_list=[val for subl in std_intraburst_freq_per_cell for val in subl]
  
	intraindices_full=nonzero(intraburst_freq)
	intraindices=sorted(set(intraindices_full[1]))
	intraburst_freq_nonzero=intraburst_freq[:,intraindices] #[1]]
	intraburst_freq_nonzero_no_trans=intraburst_freq_nonzero[intrabin_ind_no_trans:-1][:]
	num_cells_w_intraBF=len(intraindices) #[1]) #len(fullburst_freq_nonzero_no_trans)
	
	avg_interBF_per_bin_list=[val for subl in avg_interburst_freq_per_bin for val in subl]
	std_interBF_per_bin_list=[val for subl in std_interburst_freq_per_bin for val in subl]
	avg_interBF_per_cell_list=[val for subl in avg_interburst_freq_per_cell for val in subl]
	std_interBF_per_cell_list=[val for subl in std_interburst_freq_per_cell for val in subl]
  
	interindices_full=nonzero(interburst_freq)
	interindices=sorted(set(interindices_full[1]))
	interburst_freq_nonzero=interburst_freq[:,interindices] #[1]]
	interburst_freq_nonzero_no_trans=interburst_freq_nonzero[interbin_ind_no_trans:-1][:]
	num_cells_w_interBF=len(interindices) #[1]) #len(fullburst_freq_nonzero_no_trans)

   
	return [indices,num_cells_w_FBF,fullburst_freq_nonzero,fullburst_freq_nonzero_no_trans, 
			num_cells_spk_per_bin,num_cells_spk_per_bin_no_trans,avg_num_cells_spk,
			std_num_cells_spk,num_spks_per_cell_per_bin,num_spks_per_cell_per_bin_no_trans,
			avg_num_spks_per_cell,std_num_spks_per_cell,tot_num_spks_per_bin,
			tot_num_spks_per_bin_no_trans,avg_tot_num_spks,std_tot_num_spks,avg_fullburst_freq_per_bin_tot,
			std_fullburst_freq_per_bin_tot,avg_fullburst_freq_per_cell_tot,std_fullburst_freq_per_cell_tot,
			avg_fullburst_freq_per_bin,std_fullburst_freq_per_bin,
			avg_FBF_per_bin_list,std_FBF_per_bin_list,avg_FBF_per_cell_list,std_FBF_per_cell_list,
			intraindices,num_cells_w_intraBF,intraburst_freq_nonzero,intraburst_freq_nonzero_no_trans, 
			num_cells_spk_per_intrabin,num_cells_spk_per_intrabin_no_trans,avg_num_cells_spk_intra,
			std_num_cells_spk_intra,num_spks_per_cell_per_intrabin,num_spks_per_cell_per_intrabin_no_trans,
			avg_num_spks_per_cell_intra,std_num_spks_per_cell_intra,tot_num_spks_per_intrabin,
			tot_num_spks_per_intrabin_no_trans,avg_tot_num_spks_intra,std_tot_num_spks_intra,avg_intraburst_freq_per_bin_tot,
			std_intraburst_freq_per_bin_tot,avg_intraburst_freq_per_cell_tot,std_intraburst_freq_per_cell_tot,
			avg_intraburst_freq_per_bin,std_intraburst_freq_per_bin,
			avg_intraBF_per_bin_list,std_intraBF_per_bin_list,avg_intraBF_per_cell_list,std_intraBF_per_cell_list,
			interindices,num_cells_w_interBF,interburst_freq_nonzero,interburst_freq_nonzero_no_trans, 
			num_cells_spk_per_interbin,num_cells_spk_per_interbin_no_trans,avg_num_cells_spk_inter,
			std_num_cells_spk_inter,num_spks_per_cell_per_interbin,num_spks_per_cell_per_interbin_no_trans,
			avg_num_spks_per_cell_inter,std_num_spks_per_cell_inter,tot_num_spks_per_interbin,
			tot_num_spks_per_interbin_no_trans,avg_tot_num_spks_inter,std_tot_num_spks_inter,avg_interburst_freq_per_bin_tot,
			std_interburst_freq_per_bin_tot,avg_interburst_freq_per_cell_tot,std_interburst_freq_per_cell_tot,
			avg_interburst_freq_per_bin,std_interburst_freq_per_bin,
			avg_interBF_per_bin_list,std_interBF_per_bin_list,avg_interBF_per_cell_list,std_interBF_per_cell_list]
	
################################################

def make_stats_plots(num_popln,N,ctr_of_bin_ms,ctr_of_intrabin_ms,ctr_of_interbin_ms,transient,
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
					 celltype,filename,foldername):
	
	
	for i in xrange(num_popln):
		
		##Number of cells spiking per bin
		figure(1)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_bin_ms,num_cells_spk_per_bin[i])
		xlim([transient,ctr_of_bin_ms[-2]])  #skip last bin because can be incomplete (goes from second last bin to end, so info may not be helpful)
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_num_cells_spk[i], std_num_cells_spk[i]))
		
		##Number of cells spiking per intrabin
		figure(2)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_intrabin_ms,num_cells_spk_per_intrabin[i])
		xlim([transient,ctr_of_intrabin_ms[-1]])   #go to last bin since last one counted for intra and inter is complete bin  
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_num_cells_spk_intra[i], std_num_cells_spk_intra[i]))
		
		##Number of cells spiking per interbin
		figure(3)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_interbin_ms,num_cells_spk_per_interbin[i])
		xlim([transient,ctr_of_interbin_ms[-1]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_num_cells_spk_inter[i], std_num_cells_spk_inter[i]))
		
		###############        
		##Proportion of cells spiking per bin
		figure(4)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_bin_ms,num_cells_spk_per_bin[i]/float(N[i]))
		xlim([transient,ctr_of_bin_ms[-2]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_cells_spk[i]/float(N[i]), std_num_cells_spk[i]/float(N[i])))
	
		##Proportion of cells spiking per intrabin
		figure(5)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_intrabin_ms,num_cells_spk_per_intrabin[i]/float(N[i]))
		xlim([transient,ctr_of_intrabin_ms[-1]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_cells_spk_intra[i]/float(N[i]), std_num_cells_spk_intra[i]/float(N[i])))
		
		##Proportion of cells spiking per interbin
		figure(6)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_interbin_ms,num_cells_spk_per_interbin[i]/float(N[i]))
		xlim([transient,ctr_of_interbin_ms[-1]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_cells_spk_inter[i]/float(N[i]), std_num_cells_spk_inter[i]/float(N[i])))
		
		############
		###Number of spikes per cell per bin
		figure(7)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_bin_ms,num_spks_per_cell_per_bin[i])
		xlim([transient,ctr_of_bin_ms[-2]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_spks_per_cell[i],std_num_spks_per_cell[i]))
		
		###Number of spikes per cell per intrabin
		figure(8)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_intrabin_ms,num_spks_per_cell_per_intrabin[i])
		xlim([transient,ctr_of_intrabin_ms[-1]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_spks_per_cell_intra[i],std_num_spks_per_cell_intra[i]))
		
		###Number of spikes per cell per interbin
		figure(9)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_interbin_ms,num_spks_per_cell_per_interbin[i])
		xlim([transient,ctr_of_interbin_ms[-1]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_spks_per_cell_inter[i],std_num_spks_per_cell_inter[i]))
		
		########################
		### Total number of spikes per bin
		figure(10)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_bin_ms,tot_num_spks_per_bin[i])
		xlim([transient,ctr_of_bin_ms[-2]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_tot_num_spks[i],std_tot_num_spks[i]))
		
		### Total number of spikes per intrabin
		figure(11)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_intrabin_ms,tot_num_spks_per_intrabin[i])
		xlim([transient,ctr_of_intrabin_ms[-1]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_tot_num_spks_intra[i],std_tot_num_spks_intra[i]))
		
		### Total number of spikes per interbin
		figure(12)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		plot(ctr_of_interbin_ms,tot_num_spks_per_interbin[i])
		xlim([transient,ctr_of_interbin_ms[-1]])
		ylabel('%s'%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_tot_num_spks_inter[i],std_tot_num_spks_inter[i]))
		
		############
		## Fullburst Freq per bin
		figure(13)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		avg_plus_std_FBF_per_bin_list=[val for subl in (avg_fullburst_freq_per_bin[i]+std_fullburst_freq_per_bin[i]) for val in subl]
		avg_minus_std_FBF_per_bin_list=[val for subl in (avg_fullburst_freq_per_bin[i]-std_fullburst_freq_per_bin[i]) for val in subl]
		plot(ctr_of_bin_ms,avg_FBF_per_bin_list[i])
		plot(ctr_of_bin_ms,avg_plus_std_FBF_per_bin_list,'--r')
		plot(ctr_of_bin_ms,avg_minus_std_FBF_per_bin_list,'--r')
		xlim([transient,ctr_of_bin_ms[-2]])
		ylabel("%s"%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_fullburst_freq_per_bin_tot[i],std_fullburst_freq_per_bin_tot[i]))
		
		## Intraburst Freq per bin
		figure(14)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		avg_plus_std_intraBF_per_bin_list=[val for subl in (avg_intraburst_freq_per_bin[i]+std_intraburst_freq_per_bin[i]) for val in subl]
		avg_minus_std_intraBF_per_bin_list=[val for subl in (avg_intraburst_freq_per_bin[i]-std_intraburst_freq_per_bin[i]) for val in subl]
		plot(ctr_of_intrabin_ms,avg_intraBF_per_bin_list[i])
		plot(ctr_of_intrabin_ms,avg_plus_std_intraBF_per_bin_list,'--r')
		plot(ctr_of_intrabin_ms,avg_minus_std_intraBF_per_bin_list,'--r')
		xlim([transient,ctr_of_intrabin_ms[-1]])
		ylabel("%s"%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_intraburst_freq_per_bin_tot[i],std_intraburst_freq_per_bin_tot[i]))
		
		## Intraburst Freq per bin
		figure(15)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		avg_plus_std_interBF_per_bin_list=[val for subl in (avg_interburst_freq_per_bin[i]+std_interburst_freq_per_bin[i]) for val in subl]
		avg_minus_std_interBF_per_bin_list=[val for subl in (avg_interburst_freq_per_bin[i]-std_interburst_freq_per_bin[i]) for val in subl]
		plot(ctr_of_interbin_ms,avg_interBF_per_bin_list[i])
		plot(ctr_of_interbin_ms,avg_plus_std_interBF_per_bin_list,'--r')
		plot(ctr_of_interbin_ms,avg_minus_std_interBF_per_bin_list,'--r')
		xlim([transient,ctr_of_interbin_ms[-1]])
		ylabel("%s"%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_interburst_freq_per_bin_tot[i],std_interburst_freq_per_bin_tot[i]))
		
		###################
		## Fullburst Freq per cell
		figure(16)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		errorbar(xrange(N[i]),avg_FBF_per_cell_list[i],std_FBF_per_cell_list[i],linestyle='None',marker='^')
		xlim([0,100])
		ylabel("%s"%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_fullburst_freq_per_cell_tot[i],std_fullburst_freq_per_cell_tot[i]))
		
		## Intraburst Freq per cell
		figure(17)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		errorbar(xrange(N[i]),avg_intraBF_per_cell_list[i],std_intraBF_per_cell_list[i],linestyle='None',marker='^')
		xlim([0,100])
		ylabel("%s"%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_intraburst_freq_per_cell_tot[i],std_intraburst_freq_per_cell_tot[i]))
		
		## Interburst Freq per cell
		figure(18)
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(num_popln/2,2,i+1)
		errorbar(xrange(N[i]),avg_interBF_per_cell_list[i],std_interBF_per_cell_list[i],linestyle='None',marker='^')
		xlim([0,100])
		ylabel("%s"%(celltype[i]))
		title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_interburst_freq_per_cell_tot[i],std_interburst_freq_per_cell_tot[i]))
	
	############
	#### add in axes and save
	figure(1)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Number of Cells Spiking Per Full Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-num_cells_spk_per_bin.png")
	close()
	
	figure(2)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Number of Cells Spiking Per Intraburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-num_cells_spk_per_intrabin.png")
	close()
	
	figure(3)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Number of Cells Spiking Per Interburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-num_cells_spk_per_interbin.png")
	close()
	######
	
	figure(4)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Proportion of Cells Spiking Per Full Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-prop_cells_spk_per_bin.png")
	close()
	
	figure(5)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Proportion of Cells Spiking Per Intraburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-prop_cells_spk_per_intrabin.png")
	close()
	
	figure(6)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Proportion of Cells Spiking Per Interburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-prop_cells_spk_per_interbin.png")
	close()
	
	###############3
	figure(7)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Average Number of Spikes Per Cell Per Full Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-num_spk_per_cell_per_bin.png")
	close()
	
	figure(8)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Average Number of Spikes Per Cell Per Intraburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-num_spk_per_cell_per_intrabin.png")
	close()
	
	figure(9)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Average Number of Spikes Per Cell Per Interburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-num_spk_per_cell_per_interbin.png")
	close()
	
	##################
	
	figure(10)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Total Number of Spikes Per Full Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-tot_num_spks_per_bin.png")
	close()
	
	figure(11)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Total Number of Spikes Per Intraburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-tot_num_spks_per_intrabin.png")
	close()
	
	figure(12)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Total Number of Spikes Per Interburst Bin")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-tot_num_spks_per_interbin.png")
	close()
	
	#################
	figure(13)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Frequency Per Full Bin (Hz)")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-FBF_per_bin.png")
	close()
	
	figure(14)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Intraburst Frequency Per burst Bin (Hz)")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-intraBF_per_bin.png")
	close()
	
	figure(15)
	if num_popln<4:
		xlabel("Time (ms)")
	else:
		xlabel("Time (ms)")
		subplot(num_popln/2,2,3)
		xlabel("Time (ms)")
	suptitle("Interburst Frequency Per Bin (Hz)")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-interBF_per_bin.png")
	close()
	
	##############
	figure(16)
	if num_popln<4:
		xlabel("Cell Number")
	else:
		xlabel("Cell Number")
		subplot(num_popln/2,2,3)
		xlabel("Cell Number")
	suptitle("Burst Frequency Per Cell(Hz)")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-FBF_per_cell.png")
	close()
	
	figure(17)
	if num_popln<4:
		xlabel("Cell Number")
	else:
		xlabel("Cell Number")
		subplot(num_popln/2,2,3)
		xlabel("Cell Number")
	suptitle("Intraburst Frequency Per Cell (Hz)")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-intraBF_per_cell.png")
	close()
	
	figure(18)
	if num_popln<4:
		xlabel("Cell Number")
	else:
		xlabel("Cell Number")
		subplot(num_popln/2,2,3)
		xlabel("Cell Number")
	suptitle("Interburst Frequency Per Cell (Hz)")
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"-interBF_per_cell.png")
	close()
	
	
	#Dynamic frequency calculation and plotting - made by Anton
	ctr_of_bin = array(ctr_of_bin_ms) / 1000.0

	#len = 8
	# 0 1 2 3 4 5 6 7 
	dynamic_freq = zeros(len(ctr_of_bin) - 5)
	dynamic_freq_t = zeros(len(ctr_of_bin) - 5)
	

	for index in range(2,len(ctr_of_bin) - 3):
		dynamic_freq[index - 2] = 4.0/(ctr_of_bin[index + 2] - ctr_of_bin[index - 2])
		dynamic_freq_t[index - 2] = ctr_of_bin[index]
	
	#print(dynamic_freq)
	#print(dynamic_freq_t)
	
	figure(19)
	plot(dynamic_freq_t,dynamic_freq)
	title("Average burst frequency as a function of time")
	ylabel('Hz')
	xlabel('Time (s)')
	savefig(foldername+"/"+filename+"-moving_freq.png")
	close()

#'''
#Created on 2014-05-02
#
#@author: kferguso
#'''
#from scipy.signal import filtfilt,butter
#from brian.stdunits import ms
#from numpy.ma.core import zeros,ceil,logical_and, where, reshape, floor, argmax, mean, std, exp, nonzero
#import bisect
#from pylab import average
#from matplotlib.pyplot import figure, plot, ylim, xlim, ylabel, xlabel, title, show, savefig, close, suptitle, errorbar, subplot,tight_layout
#from matplotlib.ticker import MultipleLocator
#
###### determine stats: the number of spikes per burst, number of cells spiking, etc.
###collect all stats, but find averages, etc by eliminating transient
#def bursting_stats(bin_s,bin_ms,binwidth_ms,numbins,N,M,transient):
#          
#    numspks=zeros((numbins,N))
#    fullburst_freq=zeros((numbins,N))
#    
#    for i in xrange(numbins):
#        step_start=bin_s[i]
#        step_end=bin_s[i+1]
#        for j in xrange(N):
#            spks=where(logical_and(M[j]>step_start, M[j]<step_end))
#            numspks[i,j]=len(M[j][spks])
#            #find intraburst freq
#            if numspks[i,j]>1:
#                spktms=M[j][spks]
#                fullburst_freq[i,j]=1/((spktms[1:]-spktms[:-1]).mean(axis=0))
#           
#    num_cells_spk_per_bin=(numspks!=0).sum(1)
#    
#    #find averages and stds after transient
#    bin_ind_no_trans=bisect.bisect(bin_ms, transient)
#    
#    num_cells_spk_per_bin_no_trans=num_cells_spk_per_bin[bin_ind_no_trans:-2]  #ignore transient and last bin
#    avg_num_cells_spk=num_cells_spk_per_bin_no_trans.mean(axis=0)
#    std_num_cells_spk=num_cells_spk_per_bin_no_trans.std(axis=0)
#    
#    num_spks_per_cell_per_bin=numspks.mean(axis=1)
#    num_spks_per_cell_per_bin_no_trans=num_spks_per_cell_per_bin[bin_ind_no_trans:-2]
#    avg_num_spks_per_cell=num_spks_per_cell_per_bin_no_trans.mean(axis=0)
#    std_num_spks_per_cell=num_spks_per_cell_per_bin_no_trans.std(axis=0)
#    
#    tot_num_spks_per_bin=numspks.sum(1)
#    tot_num_spks_per_bin_no_trans=tot_num_spks_per_bin[bin_ind_no_trans:-2]
#    avg_tot_num_spks=tot_num_spks_per_bin.mean(axis=0)
#    std_tot_num_spks=tot_num_spks_per_bin.std(axis=0)
#    
#    avg_fullburst_freq_per_bin=zeros((numbins,1))   #can do these stats on all bins and ignore trans for total average 
#    std_fullburst_freq_per_bin=zeros((numbins,1))    
#    avg_fullburst_freq_per_cell=zeros((N,1))      #must ignore trans immediately here... 
#    std_fullburst_freq_per_cell=zeros((N,1))
#    
#    FBF_ind=nonzero(fullburst_freq)
#    fullburst_freq_no_trans=fullburst_freq[bin_ind_no_trans:-2][:]
#    FBF_ind_no_trans=nonzero(fullburst_freq_no_trans)
#    uniq_ind_per_bin=sorted(set(FBF_ind[0]))    #rows = bins
#    uniq_ind_per_cell=sorted(set(FBF_ind_no_trans[1]))    #colms = cell
#    
#    for i in uniq_ind_per_bin:
#        nonzero_row_ind=[j for j, x in enumerate(FBF_ind[0]) if x==i]
#        FBF_nonzero_row=[]
#        for k in xrange(len(nonzero_row_ind)):
#            FBF_nonzero_row.append(fullburst_freq[FBF_ind[0][nonzero_row_ind[k]]][FBF_ind[1][nonzero_row_ind[k]]])
#        avg_fullburst_freq_per_bin[i]=average(FBF_nonzero_row)
#        std_fullburst_freq_per_bin[i]=std(FBF_nonzero_row)
#         
#    
#    for i in uniq_ind_per_cell:
#        nonzero_col_ind=[j for j, x in enumerate(FBF_ind_no_trans[1]) if x==i]
#        FBF_nonzero_col_no_trans=[]
#        for k in xrange(len(nonzero_col_ind)):
#            FBF_nonzero_col_no_trans.append(fullburst_freq_no_trans[FBF_ind_no_trans[0][nonzero_col_ind[k]]][FBF_ind_no_trans[1][nonzero_col_ind[k]]])
#        avg_fullburst_freq_per_cell[i]=average(FBF_nonzero_col_no_trans)
#        std_fullburst_freq_per_cell[i]=std(FBF_nonzero_col_no_trans) 
#    
#    #Find avg FBF. Only find average of bins that have an FBF (i.e. eliminate non-zero bins).  Here we have found the average of the average FBF/bin.  That is, we have represented each bin equally, rather than each set of spikes....
#    avg_fullburst_freq_per_bin_no_trans=avg_fullburst_freq_per_bin[bin_ind_no_trans:-2]   #must eliminate trans here since didn't before
#    if max(avg_fullburst_freq_per_bin_no_trans)>0:    
#        avg_fullburst_freq_per_bin_tot=(avg_fullburst_freq_per_bin_no_trans[nonzero(avg_fullburst_freq_per_bin_no_trans)]).mean(axis=0)
#        std_fullburst_freq_per_bin_tot=(avg_fullburst_freq_per_bin_no_trans[nonzero(avg_fullburst_freq_per_bin_no_trans)]).std(axis=0)
#    else:
#        avg_fullburst_freq_per_bin_tot=0
#        std_fullburst_freq_per_bin_tot=0
#    
#    if max(avg_fullburst_freq_per_cell)>0:    
#        avg_fullburst_freq_per_cell_tot=(avg_fullburst_freq_per_cell[nonzero(avg_fullburst_freq_per_cell)]).mean(axis=0)
#        std_fullburst_freq_per_cell_tot=(avg_fullburst_freq_per_cell[nonzero(avg_fullburst_freq_per_cell)]).std(axis=0)
#    else:
#        avg_fullburst_freq_per_cell_tot=0
#        std_fullburst_freq_per_cell_tot=0
#    
#    #convert FBF to list for writing and plotting
#    avg_FBF_per_bin_list=[val for subl in avg_fullburst_freq_per_bin for val in subl]
#    std_FBF_per_bin_list=[val for subl in std_fullburst_freq_per_bin for val in subl]
#    avg_FBF_per_cell_list=[val for subl in avg_fullburst_freq_per_cell for val in subl]
#    std_FBF_per_cell_list=[val for subl in std_fullburst_freq_per_cell for val in subl]
#  
#    indices_full=nonzero(fullburst_freq)
#    indices=sorted(set(indices_full[1]))
#    fullburst_freq_nonzero=fullburst_freq[:,indices] #[1]]
#    fullburst_freq_nonzero_no_trans=fullburst_freq_nonzero[bin_ind_no_trans:-2][:]
#    num_cells_w_FBF=len(indices) #[1]) #len(fullburst_freq_nonzero_no_trans)
#   
#    return [indices,num_cells_w_FBF,fullburst_freq_nonzero,fullburst_freq_nonzero_no_trans, 
#            bin_ind_no_trans, num_cells_spk_per_bin,
#            num_cells_spk_per_bin_no_trans,avg_num_cells_spk,
#            std_num_cells_spk,num_spks_per_cell_per_bin,num_spks_per_cell_per_bin_no_trans,
#            avg_num_spks_per_cell,std_num_spks_per_cell,tot_num_spks_per_bin,
#            tot_num_spks_per_bin_no_trans,avg_tot_num_spks,std_tot_num_spks,avg_fullburst_freq_per_bin_tot,
#            std_fullburst_freq_per_bin_tot,avg_fullburst_freq_per_cell_tot,std_fullburst_freq_per_cell_tot,
#            avg_fullburst_freq_per_bin,std_fullburst_freq_per_bin,
#            avg_FBF_per_bin_list,std_FBF_per_bin_list,avg_FBF_per_cell_list,std_FBF_per_cell_list]
#    
#################################################
#
#def make_stats_plots(num_popln,N,ctr_of_bin_ms,transient,
#                     num_cells_spk_per_bin,avg_num_cells_spk,std_num_cells_spk,
#                     num_spks_per_cell_per_bin,avg_num_spks_per_cell,std_num_spks_per_cell,
#                     tot_num_spks_per_bin,avg_tot_num_spks,std_tot_num_spks,
#                     avg_fullburst_freq_per_bin,std_fullburst_freq_per_bin,
#                     avg_FBF_per_bin_list,
#                     avg_fullburst_freq_per_bin_tot,std_fullburst_freq_per_bin_tot,
#                     avg_FBF_per_cell_list,std_FBF_per_cell_list,
#                     avg_fullburst_freq_per_cell_tot,std_fullburst_freq_per_cell_tot,
#                     celltype,filename,foldername):
#    
#    
#    for i in xrange(num_popln):
#        
#        ##Number of cells spiking per bin
#    
#        figure(1)
#        if num_popln<4:
#            subplot(num_popln,1,i+1)
#        else:
#            subplot(num_popln/2,2,i+1)
#        plot(ctr_of_bin_ms,num_cells_spk_per_bin[i])
#        xlim([transient,ctr_of_bin_ms[-2]])
#        ylabel('%s'%(celltype[i]))
#        title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_num_cells_spk[i], std_num_cells_spk[i]))
#                
#        ##Proportion of cells spiking per bin
#        figure(2)
#        if num_popln<4:
#            subplot(num_popln,1,i+1)
#        else:
#            subplot(num_popln/2,2,i+1)
#        plot(ctr_of_bin_ms,num_cells_spk_per_bin[i]/float(N[i]))
#        xlim([transient,ctr_of_bin_ms[-2]])
#        ylabel('%s'%(celltype[i]))
#        title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_cells_spk[i]/float(N[i]), std_num_cells_spk[i]/float(N[i])))
#                        
#        ###Number of spikes per cell per bin
#        figure(3)
#        if num_popln<4:
#            subplot(num_popln,1,i+1)
#        else:
#            subplot(num_popln/2,2,i+1)
#        plot(ctr_of_bin_ms,num_spks_per_cell_per_bin[i])
#        xlim([transient,ctr_of_bin_ms[-2]])
#        ylabel('%s'%(celltype[i]))
#        title("%s, Avg=%0.2f, SD=%0.2f"%(celltype[i],avg_num_spks_per_cell[i],std_num_spks_per_cell[i]))
#                
#        ### Total number of spikes per bin
#        figure(4)
#        if num_popln<4:
#            subplot(num_popln,1,i+1)
#        else:
#            subplot(num_popln/2,2,i+1)
#        plot(ctr_of_bin_ms,tot_num_spks_per_bin[i])
#        xlim([transient,ctr_of_bin_ms[-2]])
#        ylabel('%s'%(celltype[i]))
#        title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_tot_num_spks[i],std_tot_num_spks[i]))
#        
#        ## Intraburst Freq per bin
#        figure(5)
#        if num_popln<4:
#            subplot(num_popln,1,i+1)
#        else:
#            subplot(num_popln/2,2,i+1)
#        avg_plus_std_FBF_per_bin_list=[val for subl in (avg_fullburst_freq_per_bin[i]+std_fullburst_freq_per_bin[i]) for val in subl]
#        avg_minus_std_FBF_per_bin_list=[val for subl in (avg_fullburst_freq_per_bin[i]-std_fullburst_freq_per_bin[i]) for val in subl]
#        plot(ctr_of_bin_ms,avg_FBF_per_bin_list[i])
#        plot(ctr_of_bin_ms,avg_plus_std_FBF_per_bin_list,'--r')
#        plot(ctr_of_bin_ms,avg_minus_std_FBF_per_bin_list,'--r')
#        xlim([transient,ctr_of_bin_ms[-2]])
#        ylabel("%s"%(celltype[i]))
#        title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_fullburst_freq_per_bin_tot[i],std_fullburst_freq_per_bin_tot[i]))
#                
#        ## Intraburst Freq per cell
#        figure(6)
#        if num_popln<4:
#            subplot(num_popln,1,i+1)
#        else:
#            subplot(num_popln/2,2,i+1)
#        errorbar(xrange(N[i]),avg_FBF_per_cell_list[i],std_FBF_per_cell_list[i],linestyle='None',marker='^')
#        xlim([0,100])
#        ylabel("%s"%(celltype[i]))
#        title("%s, Avg=%0.1f, SD=%0.1f"%(celltype[i],avg_fullburst_freq_per_cell_tot[i],std_fullburst_freq_per_cell_tot[i]))
#                
#    figure(1)
#    if num_popln<4:
#        xlabel("Time (ms)")
#    else:
#        xlabel("Time (ms)")
#        subplot(num_popln/2,2,3)
#        xlabel("Time (ms)")
#    suptitle("Number of Cells Spiking Per Bin")
#    tight_layout(pad=2.5)
#    savefig(foldername+"/"+filename+"-num_cells_spk_per_bin.png")
#    close()
#    
#    figure(2)
#    if num_popln<4:
#        xlabel("Time (ms)")
#    else:
#        xlabel("Time (ms)")
#        subplot(num_popln/2,2,3)
#        xlabel("Time (ms)")
#    suptitle("Proportion of Cells Spiking Per Bin")
#    tight_layout(pad=2.5)
#    savefig(foldername+"/"+filename+"-prop_cells_spk_per_bin.png")
#    close()
#    
#    figure(3)
#    if num_popln<4:
#        xlabel("Time (ms)")
#    else:
#        xlabel("Time (ms)")
#        subplot(num_popln/2,2,3)
#        xlabel("Time (ms)")
#    suptitle("Average Number of Spikes Per Cell Per Bin")
#    tight_layout(pad=2.5)
#    savefig(foldername+"/"+filename+"-num_spk_per_cell_per_bin.png")
#    close()
#    
#    figure(4)
#    if num_popln<4:
#        xlabel("Time (ms)")
#    else:
#        xlabel("Time (ms)")
#        subplot(num_popln/2,2,3)
#        xlabel("Time (ms)")
#    suptitle("Total Number of Spikes Per Bin")
#    tight_layout(pad=2.5)
#    savefig(foldername+"/"+filename+"-tot_num_spks_per_bin.png")
#    close()
#    
#    figure(5)
#    if num_popln<4:
#        xlabel("Time (ms)")
#    else:
#        xlabel("Time (ms)")
#        subplot(num_popln/2,2,3)
#        xlabel("Time (ms)")
#    suptitle("Intraburst Frequency Per Bin (Hz)")
#    tight_layout(pad=2.5)
#    savefig(foldername+"/"+filename+"-FBF_per_bin.png")
#    close()
#    
#    figure(6)
#    if num_popln<4:
#        xlabel("Cell Number")
#    else:
#        xlabel("Cell Number")
#        subplot(num_popln/2,2,3)
#        xlabel("Cell Number")
#    suptitle("Intraburst Frequency Per Cell (Hz)")
#    tight_layout(pad=2.5)
#    savefig(foldername+"/"+filename+"-FBF_per_cell.png")
#    close()
