
'''
Created on 2014-05-02

@author: kferguso
Edited by Anton Lunyov on 1 Aug, 2017 to work with Brian2
'''
from brian2 import *

from scipy.signal import filtfilt,butter
#from brian2 import *#from brian.stdunits import ms
from numpy.ma.core import zeros,ceil,logical_and, where, reshape, floor, argmax, mean, std, exp, arange
from matplotlib.pyplot import figure, plot, ylim, xlim, ylabel, xlabel, title, show, savefig, close, subplot,suptitle, tight_layout
from scipy.optimize import curve_fit
import bisect


#determine bin width based on spike distribution of PYR cells - for stats
def find_bursts(duration,dt,transient,N,M_t,M_i,max_freq):
	base=2 #round lgbinwidth to nearest 2 so will always divide into durations
	expnum=2.0264*exp(-0.2656*max_freq+2.9288)+5.7907
	lgbinwidth=(int(base*round((-max_freq+33)/base)))*ms   #23-good for higher freq stuff 
	#lgbinwidth=(int(base*round((expnum)/base)))/1000   #use exptl based on some fit of choice binwidths
	#lgbinwidth=10*ms
	
	numlgbins=int(ceil(duration/lgbinwidth))
	#totspkhist=zeros((numlgbins,1))
	totspkhist=zeros(numlgbins)
	#totspkdist_smooth=zeros((numlgbins,1))
	skiptime=transient*ms
	skipbin=int(ceil(skiptime/lgbinwidth))   
	
	inc_past_thresh=[]
	dec_past_thresh=[]
	
	#Create histogram given the bins calculated
	for i in xrange(numlgbins):
		step_start=(i)*lgbinwidth
		step_end=(i+1)*lgbinwidth
		totspkhist[i] = len(M_i[logical_and(M_t>step_start, M_t<step_end)])
	
	###smooth plot first so thresholds work better
	#totspkhist_1D=reshape(totspkhist,len(totspkhist))  #first just reshape so single row not single colm
	#b,a=butter(3,0.4,'low')
	#totspkhist_smooth=filtfilt(b,a,totspkhist_1D)
		
	
	#totspkhist_smooth=reshape(totspkhist,len(totspkhist))  #here we took out the actual smoothing and left it as raw distn. here just reshape so single row not single colm
	totspkdist_smooth=totspkhist/max(totspkhist[skipbin:])  #create distn based on hist, but skip first skiptime to cut out transient excessive spiking
	
#    ####### FOR MOVING THRESHOLD #################
	## find points where increases and decreases over some threshold
	dist_thresh=[]
	thresh_plot=[]

	mul_fac=0.35
	switch=0  #keeps track of whether inc or dec last 
	elim_noise=1/(max_freq*2.5*Hz)
	#For line 95, somehow not required in previous version?
	#elim_noise_units = 1/(max_freq*Hz*2.5)
	
	thresh_time=5/(max_freq)#capture 5 cycles
	thresh_ind=int(floor((thresh_time/lgbinwidth)/2)) 	#the number of indices on each side of the window
	
	#dist_thresh moves with window capturing approx 5 cycles (need special cases for borders) Find where increases and decreases past threshold (as long as a certain distance apart, based on "elim_noise" which is based on avg freq of bursts
	dist_thresh.append(totspkdist_smooth[skipbin:skipbin+thresh_ind].mean(0)+mul_fac*totspkdist_smooth[skipbin:skipbin+thresh_ind].std(0))
	
	for i in xrange(1,numlgbins):
		step_start=(i)*lgbinwidth
		step_end=(i+1)*lgbinwidth
		
		#moving threshold
		if i>(skipbin+thresh_ind) and (i+thresh_ind)<len(totspkdist_smooth):
			#print(totspkdist_smooth[i-thresh_ind:i+thresh_ind])
			dist_thresh.append(totspkdist_smooth[i-thresh_ind:i+thresh_ind].mean(0)+mul_fac*totspkdist_smooth[i-thresh_ind:i+thresh_ind].std(0))
		elif (i+thresh_ind)>=len(totspkdist_smooth):
			dist_thresh.append(totspkdist_smooth[-thresh_ind:].mean(0)+mul_fac*totspkdist_smooth[-thresh_ind:].std(0))
		else:
			dist_thresh.append(totspkdist_smooth[skipbin:skipbin+thresh_ind].mean(0)+mul_fac*totspkdist_smooth[skipbin:skipbin+thresh_ind].std(0))
		
		
		if (totspkdist_smooth[i-1]<dist_thresh[i]) and (totspkdist_smooth[i]>=dist_thresh[i]):
			#inc_past_thresh.append(step_start-0.5*lgbinwidth)   
			if (inc_past_thresh):    #there has already been at least one inc, 
				if (abs(inc_past_thresh[-1]-(step_start-0.5*lgbinwidth))>elim_noise) and switch==0:   #must be at least x ms apart (yHz), and it was dec last..
					inc_past_thresh.append(step_start-0.5*lgbinwidth)  #take lower point (therefore first) when increasing. Need to -0.5binwidth to adjust for shift between index of bin width and index of bin distn
					#print (['incr=%f'%inc_past_thresh[-1]])
					thresh_plot.append(dist_thresh[i])    
					switch=1            
			else:
				inc_past_thresh.append(step_start-0.5*lgbinwidth)  #take lower point (therefore first) when increasing. Need to -0.5binwidth to adjust for shift between index of bin width and index of bin distn
				thresh_plot.append(dist_thresh[i])
				switch=1   #keeps track of that it was inc. last
		elif (totspkdist_smooth[i-1]>=dist_thresh[i]) and (totspkdist_smooth[i]<dist_thresh[i]):
		   # dec_past_thresh.append(step_end-0.5*lgbinwidth)  #take lower point (therefore second) when decreasing
			if (dec_past_thresh):    #there has already been at least one dec
				if (abs(dec_past_thresh[-1]-(step_end-0.5*lgbinwidth))>elim_noise) and switch==1:    #must be at least x ms apart (y Hz), and it was inc last
					dec_past_thresh.append(step_end-0.5*lgbinwidth)  #take lower point (therefore second) when decreasing
					#print (['decr=%f'%dec_past_thresh[-1]])
					switch=0
			else:
				dec_past_thresh.append(step_end-0.5*lgbinwidth)  #take lower point (therefore second) when decreasing
				switch=0    #keeps track of that it was dec last
	
	if totspkdist_smooth[0]<dist_thresh[0]:   #if you are starting below thresh, then pop first inc.  otherwise, don't (since will decrease first)
		if inc_past_thresh:  #if list is not empty
			inc_past_thresh.pop(0)
#  

#####################################################################
#    
	######### TO DEFINE A STATIC THRESHOLD AND FIND CROSSING POINTS
	
#    dist_thresh=0.15 #static threshold  
#    switch=0  #keeps track of whether inc or dec last
#    overall_freq=3.6 #0.9
#    elim_noise=1/(overall_freq*5)#2.5)
#    
#    
#    for i in xrange(1,numlgbins):
#        step_start=(i)*lgbinwidth
#        step_end=(i+1)*lgbinwidth
#     
#        if (totspkdist_smooth[i-1]<dist_thresh) and (totspkdist_smooth[i]>=dist_thresh):   #if cross threshold (increasing)
#            if (inc_past_thresh):    #there has already been at least one inc, 
#                if (abs(dec_past_thresh[-1]-(step_start-0.5*lgbinwidth))>elim_noise) and switch==0:   #must be at least x ms apart (yHz) from the previous dec, and it was dec last..
#                    inc_past_thresh.append(step_start-0.5*lgbinwidth)  #take lower point (therefore first) when increasing. Need to -0.5binwidth to adjust for shift between index of bin width and index of bin distn
#                    #print (['incr=%f'%inc_past_thresh[-1]])     #-0.5*lgbinwidth
#                    switch=1            
#            else:
#                inc_past_thresh.append(step_start-0.5*lgbinwidth)  #take lower point (therefore first) when increasing. Need to -0.5binwidth to adjust for shift between index of bin width and index of bin distn
#                switch=1   #keeps track of that it was inc. last
#        elif (totspkdist_smooth[i-1]>=dist_thresh) and (totspkdist_smooth[i]<dist_thresh):
#            if (dec_past_thresh):    #there has already been at least one dec
#                if (abs(inc_past_thresh[-1]-(step_end-0.5*lgbinwidth))>elim_noise) and switch==1:    #must be at least x ms apart (y Hz) from the previous incr, and it was inc last
#                    dec_past_thresh.append(step_end-0.5*lgbinwidth)  #take lower point (therefore second) when decreasing
#                    #print (['decr=%f'%dec_past_thresh[-1]])
#                    switch=0
#            else:
#                dec_past_thresh.append(step_end-0.5*lgbinwidth)  #take lower point (therefore second) when decreasing
#                switch=0    #keeps track of that it was dec last
#    
#    
#    if totspkdist_smooth[0]<dist_thresh:   #if you are starting below thresh, then pop first inc.  otherwise, don't (since will decrease first)
#        if inc_past_thresh:  #if list is not empty
#            inc_past_thresh.pop(0)
	
	
################################################################
###############################################################
	
	######## DEFINE INTER AND INTRA BURSTS ########
	
	#since always start with dec, intraburst=time points from 1st inc:2nd dec, from 2nd inc:3rd dec, etc.  
	#interburst=time points from 1st dec:1st inc, from 2nd dec:2nd inc, etc. 
	
	intraburst_time_ms_compound_list=[]
	interburst_time_ms_compound_list=[]
	intraburst_bins=[]  #in seconds
	interburst_bins=[]
	
	#print(inc_past_thresh)
	if len(inc_past_thresh)<len(dec_past_thresh):   #if you end on a decrease
		for i in xrange(len(inc_past_thresh)):
			intraburst_time_ms_compound_list.append(arange(inc_past_thresh[i]/ms,dec_past_thresh[i+1]/ms,1))  #10 is timestep
			interburst_time_ms_compound_list.append(arange((dec_past_thresh[i]+dt)/ms,(inc_past_thresh[i]-dt)/ms,1))    #10 is timestep
			intraburst_bins.append(inc_past_thresh[i])
			intraburst_bins.append(dec_past_thresh[i+1])
			interburst_bins.append(dec_past_thresh[i])
			interburst_bins.append(inc_past_thresh[i])
	else: 											#if you end on an increase
		for i in xrange(len(inc_past_thresh)-1):   
			intraburst_time_ms_compound_list.append(arange(inc_past_thresh[i]/ms,dec_past_thresh[i+1]/ms,1))  #10 is timestep
			interburst_time_ms_compound_list.append(arange((dec_past_thresh[i]+dt)/ms,(inc_past_thresh[i]-dt)/ms,1))    #10 is timestep
			intraburst_bins.append(inc_past_thresh[i])
			intraburst_bins.append(dec_past_thresh[i+1])
			interburst_bins.append(dec_past_thresh[i]+dt)
			interburst_bins.append(inc_past_thresh[i]-dt)
		if dec_past_thresh and inc_past_thresh:   #if neither dec_past_thresh nor inc_past_thresh is empty
			interburst_bins.append(dec_past_thresh[-1]+dt)   #will have one more inter than intra
			interburst_bins.append(inc_past_thresh[-1]+dt)
	
	interburst_bins = interburst_bins / second
	intraburst_bins = intraburst_bins / second
	
	intraburst_time_ms=[num for elem in intraburst_time_ms_compound_list for num in elem]   #flatten list
	interburst_time_ms=[num for elem in interburst_time_ms_compound_list for num in elem]   #flatten list
	
	num_intraburst_bins=len(intraburst_bins)/2    #/2 since have both start and end points for each bin
	num_interburst_bins=len(interburst_bins)/2

	
	intraburst_bins_ms=[x * 1000 for x in intraburst_bins]
	interburst_bins_ms=[x * 1000 for x in interburst_bins]
	
	###################################### 
	#bin_s=[((inc_past_thresh-dec_past_thresh)/2+dec_past_thresh) for inc_past_thresh, dec_past_thresh in zip(inc_past_thresh,dec_past_thresh)]
	bin_s=[((x-y)/2+y) for x,y in zip(inc_past_thresh,dec_past_thresh)] / second
	
	
	binpt_ind=[int(floor(x/lgbinwidth)) for x in bin_s]
	
	
	########## FIND PEAK TO TROUGH AND SAVE VALUES  ###################
	########## CATEGORIZE BURSTING BASED ON PEAK TO TROUGH VALUES ###################
	########## DISCARD BINPTS IF PEAK TO TROUGH IS TOO SMALL ###################
	
	peaks=[]
	trough=[]
	peak_to_trough_diff=[]
	min_burst_size=0.2 #defines a burst as 0.2 or larger.
	
	for i in xrange(len(binpt_ind)-1):
		peaks.append(max(totspkdist_smooth[binpt_ind[i]:binpt_ind[i+1]]))
		trough.append(min(totspkdist_smooth[binpt_ind[i]:binpt_ind[i+1]]))
	
	peak_to_trough_diff=[max_dist-min_dist for max_dist,min_dist in zip(peaks,trough)]
	
	#to delete all bins following any <min_burst_size
	first_ind_not_burst=next((x[0] for x in enumerate(peak_to_trough_diff) if x[1] < 0.2), None)
#    if first_ind_not_burst:
#        del bin_s[first_ind_not_burst+1:]   #needs +1 since bin_s has one additional value (since counts edges)
	
	#to keep track of any bins <0.2 so can ignore in stats later
	all_ind_not_burst=[x[0] for x in enumerate(peak_to_trough_diff) if x[1] < 0.2]   #defines a burst as 0.2 or larger.
	
	bin_ms=[x*1000 for x in bin_s]
	binpt_ind=[int(floor(x/lgbinwidth)) for x in bin_s]
	
	#for moving threshold only
	thresh_plot=[]
	thresh_plot=[dist_thresh[x] for x in binpt_ind]
	
	#for static threshold
	#thresh_plot=[dist_thresh]*len(bin_ms)
#    
#    
#    bin_s=[((inc_past_thresh-dec_past_thresh)/2+dec_past_thresh) for inc_past_thresh, dec_past_thresh in zip(inc_past_thresh,dec_past_thresh)]
#    bin_ms=[x*1000 for x in bin_s]
#    thresh_plot=[]
#    binpt_ind=[int(floor(x/lgbinwidth)) for x in bin_s]
#    thresh_plot=[dist_thresh[x] for x in binpt_ind]       
#        
	binpts=xrange(int(lgbinwidth*1000/2), int(numlgbins*lgbinwidth*1000), int(lgbinwidth*1000))
	totspkhist_list=totspkhist.tolist()#[val for subl in totspkhist for val in subl] 
   

	#find first index after transient to see if have enough bins to do stats
	bin_ind_no_trans=bisect.bisect(bin_ms, transient)
	intrabin_ind_no_trans=bisect.bisect(intraburst_bins, transient/1000)  #transient to seconds      
	if intrabin_ind_no_trans % 2 != 0:   #index must be even since format is ind0=start_bin, ind1=end_bin, ind2=start_bin, .... .
		intrabin_ind_no_trans += 1
	interbin_ind_no_trans=bisect.bisect(interburst_bins, transient/1000)
	if interbin_ind_no_trans % 2 != 0:
		interbin_ind_no_trans += 1
		
	
	return [bin_s,bin_ms,binpts,totspkhist,totspkdist_smooth,dist_thresh,
			totspkhist_list,thresh_plot,binpt_ind,lgbinwidth,numlgbins,
			intraburst_bins,interburst_bins,intraburst_bins_ms,interburst_bins_ms,intraburst_time_ms,interburst_time_ms,
			num_intraburst_bins,num_interburst_bins,bin_ind_no_trans,intrabin_ind_no_trans,interbin_ind_no_trans]
	


### find distns for other cell populations besides PYR
def find_distn(lgbinwidth,numlgbins,transient,N,M_t,M_i):            
	totspkhist=zeros((numlgbins,1))
	skiptime=transient*ms
	skipbin=int(ceil(skiptime/lgbinwidth)) 
	for i in xrange(numlgbins):
		step_start=(i)*lgbinwidth  #30*ms #
		step_end=(i+1)*lgbinwidth  #30*ms #
		for j in xrange(N):
			#spks=where(logical_and(M[j]>step_start, M[j]<step_end))
			#totspkhist[i]+=len(M[j][spks])
			totspkhist[i] += len(M_i[logical_and(M_t>step_start, M_t<step_end)])
	
	#totspkhist_1D=reshape(totspkhist,len(totspkhist))
	###smooth plot first so thresholds work better
	#b,a=butter(3,0.4,'low')
	#totspkhist_smooth=filtfilt(b,a,totspkhist_1D)
	totspkhist_smooth=reshape(totspkhist,len(totspkhist))   #here we took out the actual smoothing and left it as raw distn. 
	
	#create distn based on hist, but skip first skiptime to cut out transient excessive spiking
	if max(totspkhist_smooth[skipbin:])>0:
		totspkdist_smooth=totspkhist_smooth/max(totspkhist_smooth[skipbin:])
	else: 
		totspkdist_smooth=totspkhist_smooth
	totspkhist_list=[val for subl in totspkhist for val in subl]
	return [totspkhist,totspkdist_smooth,totspkhist_list]



###################################################################################

### make dist plots - over time and one where distns are overlaid/burst.  In the overlaid one, include mean and gaussian fit of mean if possible
def make_dist_plots(duration,transient,num_popln,lgbinwidth,bin_s,bin_ms,binpts,binpt_ind,totspkdist_smooth,thresh_plot,dist_thresh,intraburst_time_ms,interburst_time_ms,max_freq,max_pwr,g_string,ge_mean_str,ge_SD_str,celltype,filename,foldername): 
	
	### spike distribution over time
	for i in xrange(num_popln):
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(2,2,i+1)
		plot(binpts,totspkdist_smooth)        
		plot(bin_ms,thresh_plot,'*r')
		if i==0: #if PYR, also plot threshold w time:
			#plot(binpts,[dist_thresh]*len(binpts),'g')
			plot(intraburst_time_ms,[dist_thresh[0]-0.025]*len(intraburst_time_ms),'.r')
			plot(interburst_time_ms,[dist_thresh[0]-0.05]*len(interburst_time_ms),'.k')
		ylim([0,1])
		ylabel('%s'%(celltype[i]))
	xlabel("Time (ms)")
	#suptitle('Spike Distn, max_freq=%0.1f, max_pwr=%0.5f, \n gpyr=%s, ge=%s, geSD=%s'%(max_freq,max_pwr,g_string,ge_mean_str,ge_SD_str))
#    suptitle('Spike Distn, max_freq=%0.1f, max_pwr=%0.5f, \n gpyr=%s, Iapp=%s, IappSD=%s'%(max_freq,max_pwr,g_string,mean_string,hetero_string))
	#tight_layout(pad=4)
	#show()
	savefig(foldername+"/"+filename+"_PYR_dist.png")
	close()
	
	##make distn overlay plots and fit with gaussian if possible
	#_make_dist_overlay_plot(num_popln,lgbinwidth,bin_s,binpt_ind,totspkdist_smooth,celltype,filename,foldername)
	

	### spike distribution up close over time
	for i in xrange(num_popln):
		if num_popln<4:
			subplot(num_popln,1,i+1)
		else:
			subplot(2,2,i+1)
		plot(binpts,totspkdist_smooth)        
		plot(bin_ms,thresh_plot,'*r')
		if i==0: #if PYR, also plot threshold w time:
			#plot(binpts,[dist_thresh]*len(binpts),'g')
			plot(intraburst_time_ms,[dist_thresh[0]-0.025]*len(intraburst_time_ms),'.r')
			plot(interburst_time_ms,[dist_thresh[0]-0.05]*len(interburst_time_ms),'.k')
		ylim([0,1])
		xlim(duration/ms - transient,duration/ms)
		ylabel('%s'%(celltype[i]))
		xlabel("Time (ms)")
	#suptitle('Spike Distn, max_freq=%0.1f, max_pwr=%0.5f, \n gpyr=%s, ge=%s, geSD=%s'%(max_freq,max_pwr,g_string,ge_mean_str,ge_SD_str))
	#tight_layout(pad=4)
	#show()
	savefig(foldername+"/"+filename+"_PYR_dist_end.png")
	close()
	
	
	#######################################################

#make binwidth plot, and include endpoints to bin_s and bin_ms.  Write out binwidths. Create ctr_of_bins
def make_binwidth_plot(duration,bin_s,bin_ms,intraburst_bins_ms,interburst_bins_ms,transient,filename,foldername):    
	#bins include first and last point (0 and duration seconds or ms)
	insert(bin_s,0,0.0)
	append(bin_s,duration/second)
	append(bin_ms,duration/ms)
	insert(bin_ms,0,0.0)
	
	#find binwidth
#    bin_ms_temp=bin_ms[:-1]
#    bin_ms_shift=bin_ms[1:]
#    binwidth_ms=[bin_ms_shift-bin_ms_temp for bin_ms_shift,bin_ms_temp in zip(bin_ms_shift,bin_ms_temp)]
#    binwidth_ms_temp=binwidth_ms
	
	binwidth_ms=[x-y for x,y in zip(bin_ms[1:],bin_ms[:-1])]
	intrabinwidth_ms=[a - b for a, b in zip(intraburst_bins_ms[1::2], intraburst_bins_ms[::2])]   #all odd indices - all even indices
	interbinwidth_ms=[a - b for a, b in zip(interburst_bins_ms[1::2], interburst_bins_ms[::2])]   #all odd indices - all even indices
	
	#find binwidth avg and std
	avg_binwidth_ms=mean(binwidth_ms)
	std_binwidth_ms=std(binwidth_ms)    
	avg_intrabinwidth_ms=mean(intrabinwidth_ms)
	std_intrabinwidth_ms=std(intrabinwidth_ms)
	avg_interbinwidth_ms=mean(interbinwidth_ms)
	std_interbinwidth_ms=std(interbinwidth_ms)
	
	#write out binwidth info.  
	#Format: 
	#avg std
	#binwidth1 binwidth2 ....... binwidthn 
	f_binwidth_ms=open(foldername+"/"+filename+"_binwidth.txt","w")
	f_binwidth_ms.write(str(avg_binwidth_ms) + " ")
	f_binwidth_ms.write(str(std_binwidth_ms) + " ")
	f_binwidth_ms.write("\n")
	f_binwidth_ms.write(' '.join(map(str,binwidth_ms))) 
						 
	#write out intrabins
	f_intrabinwidth_ms=open(foldername+"/"+filename+"_intrabinwidth.txt","w")
	f_intrabinwidth_ms.write(str(avg_intrabinwidth_ms) + " ")
	f_intrabinwidth_ms.write(str(std_intrabinwidth_ms) + " ")
	f_intrabinwidth_ms.write("\n")
	f_intrabinwidth_ms.write(' '.join(map(str,intrabinwidth_ms))) 
	
	#write out interbins
	f_interbinwidth_ms=open(foldername+"/"+filename+"_interbinwidth.txt","w")
	f_interbinwidth_ms.write(str(avg_interbinwidth_ms) + " ")
	f_interbinwidth_ms.write(str(std_interbinwidth_ms) + " ")
	f_interbinwidth_ms.write("\n")
	f_interbinwidth_ms.write(' '.join(map(str,interbinwidth_ms))) 
	
	
	#find center of bin in time
	ctr_of_bin_ms=[x-0.5*y for x,y in zip(bin_ms[1:],binwidth_ms)]
	ctr_of_intrabin_ms=[x-0.5*y for x,y in zip(intraburst_bins_ms[1::2],intrabinwidth_ms)]
	ctr_of_interbin_ms=[x-0.5*y for x,y in zip(interburst_bins_ms[1::2],interbinwidth_ms)]
	
	#find number of bins
	numbins=len(bin_s)-1
	
	### Plot Binwidths
	plot(ctr_of_bin_ms,binwidth_ms)
	if len(ctr_of_bin_ms)>1:
		if ctr_of_bin_ms[-2]>(transient):
			xlim([transient,ctr_of_bin_ms[-2]])
	xlabel("Time (ms)")
	ylabel("Bin Width (ms)")
	#suptitle("Bin Width (ms)")
	title("Avg=%0.1f, SD=%0.1f"%(avg_binwidth_ms,std_binwidth_ms))
	#tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"_binwidth.png")
	close()
	
	plot(ctr_of_intrabin_ms,intrabinwidth_ms)
	if len(ctr_of_intrabin_ms)>1:
		if ctr_of_intrabin_ms[-2]>(transient):
			xlim([transient,ctr_of_intrabin_ms[-2]])
	xlabel("Time (ms)")
	ylabel("Intraburst Bin Width (ms)")
	#suptitle("Intraburst Bin Width (ms)")
	title("Avg=%0.1f, SD=%0.1f"%(avg_intrabinwidth_ms,std_intrabinwidth_ms))
	#tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"_intraburst_binwidth.png")
	close()
	
	plot(ctr_of_interbin_ms,interbinwidth_ms)
	if len(ctr_of_interbin_ms)>1:
		if ctr_of_interbin_ms[-2]>transient:
			xlim([transient,ctr_of_interbin_ms[-2]])
	xlabel("Time (ms)")
	ylabel("Interburst Bin Width (ms)")
	#title("Interburst Bin Width (ms)")
	title("Avg=%0.1f, SD=%0.1f"%(avg_interbinwidth_ms,std_interbinwidth_ms))
	tight_layout(pad=2.5)
	savefig(foldername+"/"+filename+"_interburst_binwidth.png")
	close()

	return[binwidth_ms,intrabinwidth_ms,interbinwidth_ms,ctr_of_bin_ms,ctr_of_intrabin_ms,ctr_of_interbin_ms,numbins,bin_ms,bin_s]


##################################################################################
#### spike dist overlaid for each bin
#def _make_dist_overlay_plot(num_popln,lgbinwidth,bin_s,binpt_ind,totspkdist_smooth,celltype,filename,foldername):
	#numbins=len(bin_s)
	#
	#
	#for i in xrange(num_popln):
	#    dist_of_dist=[[] for _ in xrange(numbins-1)]
	#    max_ind=[]
	#    
	#    for j in xrange(numbins-2): #skip last bin
	#        dist_of_dist[j]=totspkdist_smooth[i][binpt_ind[j]:binpt_ind[j+1]]
	#        max_ind.append(argmax(dist_of_dist[j]))
	#        max_val=dist_of_dist[j][max_ind[j]]
	#        if max_val>0:
	#            dist_of_dist[j]=[x/max_val for x in dist_of_dist[j]]
	#        binpt_range=binpt_ind[j+1]-binpt_ind[j]
	#        figure(1)
	#        if num_popln<4:
	#            subplot(num_popln,1,i+1)
	#        else:
	#            subplot(num_popln/2,2,i+1)
	#        plot(xrange(int(-max_ind[j]*lgbinwidth*1000),int((binpt_range-max_ind[j])*lgbinwidth*1000),int(lgbinwidth*1000)),dist_of_dist[j],color = '0.75')
	#        if j==0:
	#            min_start_ind=max_ind[j]
	#            min_end_ind=binpt_range-max_ind[j]
	#        else:
	#            min_start_ind=min(min_start_ind,max_ind[j])
	#            min_end_ind=min(min_end_ind,binpt_range-max_ind[j])
	#
	#    num_ind_for_mean=min_end_ind+min_start_ind
	#    
	#    if num_ind_for_mean>1:   #only find mean and gaussian of mean  if possible
	#        tot_dist_of_dist=[0]*num_ind_for_mean
	#    
	#        for j in xrange(numbins-2):   #skip last bin
	#            dist_temp=dist_of_dist[j][max_ind[j]-min_start_ind:max_ind[j]+min_end_ind]
	#            tot_dist_of_dist=[x+y for x,y in zip(tot_dist_of_dist,dist_temp)]
	#    
	#        mean_dist_of_dist=[x/(numbins-2) for x in tot_dist_of_dist]
	#        x_pts=xrange(int(-min_start_ind*lgbinwidth*1000),int((min_end_ind)*lgbinwidth*1000),int(lgbinwidth*1000))                       #the number of data    
	#          
	#        plot(x_pts,mean_dist_of_dist,'k')
	#        
	#        #find half-width, height, and gaussian fit - only do it at least two points on either side of max and only do for range between the max of the two endpts
	#        if argmax(mean_dist_of_dist)>1 and argmax(mean_dist_of_dist)+2<len(mean_dist_of_dist):
	#            #mean_height=max(mean_dist_of_dist)-min(mean_dist_of_dist)   #based on min point.  need to base on max of min ends
	#            higher_endpt=max(mean_dist_of_dist[0],mean_dist_of_dist[-1])
	#            mean_height=max(mean_dist_of_dist)-higher_endpt
	#            half_mean_dist_of_dist=mean_height/2 + higher_endpt #min(mean_dist_of_dist)
	#            
	#            first_ind=next(x[0] for x in enumerate(mean_dist_of_dist) if x[1]>half_mean_dist_of_dist)
	#            second_ind_temp=next(x[0] for x in enumerate(mean_dist_of_dist[first_ind+1:]) if x[1]<half_mean_dist_of_dist)
	#            second_ind=second_ind_temp+first_ind+1
	#            half_width=x_pts[second_ind]-x_pts[first_ind]
	#        
	#            title('%s, avg height=%0.1f, half-width=%0.1fms'%(celltype[i],mean_height,half_width))
	#        else:   
	#            title('%s'%(celltype[i]))
	#    ################ GAUSSIAN FIT ######################
	#        mean_dist_of_dist_shift=mean_dist_of_dist-min(mean_dist_of_dist)
	#        if max(mean_dist_of_dist_shift)>0:
	#            
	#            mean_dist_of_dist_scaled=mean_dist_of_dist_shift/max(mean_dist_of_dist_shift)
	#            
	#            #find max index and the min number of points on either side of max.  Only use this range to determine fit  
	#            max_mean_dist_of_dist_scaled = max(mean_dist_of_dist_scaled)
	#            max_pt_ind=[k for k, m in enumerate(mean_dist_of_dist_scaled) if m == max_mean_dist_of_dist_scaled]
	#            max_pt_ind=max_pt_ind[0]
	#            num_pts_right=len(mean_dist_of_dist_scaled)-(max_pt_ind+1)
	#            num_pts_left=max_pt_ind
	#            num_pts_each_side=min(num_pts_right,num_pts_left)
	#            mean_DoD=mean_dist_of_dist_scaled[max_pt_ind-num_pts_each_side:max_pt_ind+num_pts_each_side+1]
	#            mean_DoD_shift=mean_DoD-min(mean_DoD)
	#            mean_DoD_scaled=mean_DoD_shift/max(mean_DoD_shift)
	#            ##estimate mean and std for gaussian fit 
	#            x_pts_scaled=xrange(int(-num_pts_each_side*lgbinwidth*1000),int((num_pts_each_side+1)*lgbinwidth*1000),int(lgbinwidth*1000))  
	#            n = len(x_pts_scaled)   
	#            mean_of_means=mean(x_pts_scaled)
	#            sigma_of_means = std(x_pts_scaled) 
	#            
	#            def gaus(x, *p):
	#                A,mu,sigma=p
	#                return A*exp(-(x-mu)**2/(2*sigma**2))
	#            
	#            try:
	#                coeff,var_matrix = curve_fit(gaus,x_pts_scaled,mean_DoD_scaled,p0=[1,mean_of_means,sigma_of_means])   #p0 is initial guess for fitting coefficients
	#                gaus_fit=gaus(x_pts_scaled,*coeff)
	#                
	#                # get the fitting parameters, i.e. the mean and standard deviation:
	#                gaus_mean=coeff[1]
	#                gaus_std=coeff[2]
	#                
	#                figure(2)
	#                if num_popln<4:
	#                    subplot(num_popln,1,i+1)
	#                else:
	#                    subplot(num_popln/2,2,i+1)                
	#                plot(x_pts,mean_dist_of_dist_scaled,color = '0.75')
	#                plot(x_pts_scaled,mean_DoD_scaled,'k')
	#                plot(x_pts_scaled,gaus_fit,'r')
	#                title('%s, mean=%0.1f, std=%0.1f'%(celltype[i],gaus_mean,gaus_std))
	#                           
	#            except RuntimeError:
	#                print("Error - curve_fit failed")
	#    
	#    else:
	#        figure(1)
	#        if num_popln<4:
	#            subplot(num_popln,1,i+1)
	#        else:
	#            subplot(num_popln/2,2,i+1)
	#        title('%s'%(celltype[i]))
	#                
	#figure(1)
	#suptitle('Indiv (grey) and avg (black) distn')
	#tight_layout(pad=2)
	#savefig(foldername+"/"+filename+"-overlay-dist.png")
	#close()
	#figure(2)
	#suptitle('Gaussian Fit (red) of Mean Distn')
	#savefig(foldername+"/"+filename+"-gauss_fit_of_distn.png")
	#close()


