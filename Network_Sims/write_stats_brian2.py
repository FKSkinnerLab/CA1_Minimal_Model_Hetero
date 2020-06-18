'''
Created on 2014-05-03

@author: kferguso
Edited by Anton Lunyov on 3 Aug, 2017 to work with Brian2
'''

####write out stats###
def write_stats(bin_ms,binpts,num_popln,avg_num_cells_spk,std_num_cells_spk,
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
                num_interburst_bins,interburst_freq_nonzero,filename,foldername):   
       
    #Format: PYRtotavg PYRtotstd PYRbin1 PYRbin2 PYRbin3 .... 
    #        PVtotavg PVtotstd PVbin1  PVbin2  PVbin3 .....
    #        ......
    # 
    #For intraburst freq per bin, write out in separate file.  Note here write for ALL bins (incl transient and last bin), but only plot from after transient to bin(end-1)    
    #    PYR_tot_avg PYR_tot_std
    #    PYRbin1_avg PYRbin2_avg  .........   
    #    PYRbin1_std PYRbin2_std .........
    #
    #....
    #For intraburst freq per cell, write out in separate file.  Format for all values of FBF - only print out cells with at least one non-zero FBF, and first number is cell number.  
    #Note, here am not writing out transient values, but are keeping last bin value.  Plot and avg/std stats ignore trans + last bin..  
    #    num_PYR_tot num_PV_tot
    #
    #    num_PYR_tot PYRavg PYRstd
    #    num_PV_tot  PVavg  PVstd
    #    .....
    #
    #    PYRnum1  PYR1avg PYR1std bin1  bin2  ..... 
    #    PYRnum2  PYRs2avg PYR2std bin1  bin2 ....
    #    .....
    #
    #For hist/dist:
    # bin_ms1 bin_ms2 .....
    # xrange(int(lgbinwidth*1000/2), int(numlgbins*lgbinwidth*1000), int(lgbinwidth*1000))
    #
    # totPYRspkdist1 totPYRspkdist1 totPYRspkdist3 .....
    # totBCspkdist1 totBCspkdist1 totBCspkdist3 .....
    
    
    f_num_cells_spking_per_bin=open(foldername+"/"+filename+"-num_cells_spking_per_bin.txt","w")
    f_num_spks_per_cell_per_bin=open(foldername+"/"+filename+"-num_spks_per_cell_per_bin.txt","w")
    f_tot_num_spks_per_bin=open(foldername+"/"+filename+"-tot_num_spks_per_bin.txt","w")
    f_FBF_per_bin=open(foldername+"/"+filename+"-FBF_per_bin.txt","w")
    f_FBF_per_cell=open(foldername+"/"+filename+"-FBF_per_cell.txt","w")
    f_spk_hist_and_dist=open(foldername+"/"+filename+"-spk_hist_and_dist.txt","w")
    f_bins_ms=open(foldername+"/"+filename+"-bins_ms.txt","w")
    
    f_num_cells_spking_per_intrabin=open(foldername+"/"+filename+"-num_cells_spking_per_intrabin.txt","w")
    f_num_spks_per_cell_per_intrabin=open(foldername+"/"+filename+"-num_spks_per_cell_per_intrabin.txt","w")
    f_tot_num_spks_per_intrabin=open(foldername+"/"+filename+"-tot_num_spks_per_intrabin.txt","w")
    f_intraBF_per_bin=open(foldername+"/"+filename+"-intraBF_per_bin.txt","w")
    f_intraBF_per_cell=open(foldername+"/"+filename+"-intraBF_per_cell.txt","w")
    f_intraburst_bins_ms=open(foldername+"/"+filename+"-intraburst_bins_ms.txt","w")
    
    f_num_cells_spking_per_interbin=open(foldername+"/"+filename+"-num_cells_spking_per_interbin.txt","w")
    f_num_spks_per_cell_per_interbin=open(foldername+"/"+filename+"-num_spks_per_cell_per_interbin.txt","w")
    f_tot_num_spks_per_interbin=open(foldername+"/"+filename+"-tot_num_spks_per_interbin.txt","w")
    f_interBF_per_bin=open(foldername+"/"+filename+"-interBF_per_bin.txt","w")
    f_interBF_per_cell=open(foldername+"/"+filename+"-interBF_per_cell.txt","w")
    f_interburst_bins_ms=open(foldername+"/"+filename+"-interburst_bins_ms.txt","w")
       
    f_bins_ms.write(' '.join(map(str,bin_ms)))
    f_intraburst_bins_ms.write(' '.join(map(str,intraburst_bins_ms)))
    f_interburst_bins_ms.write(' '.join(map(str,interburst_bins_ms)))
    
    #for i in xrange(num_popln): 
    f_spk_hist_and_dist.write(' '.join(map(str,binpts)))
    
    f_spk_hist_and_dist.write("\n\n")
    #write out bin details for each cell. 
    for i in xrange(num_popln):
        #write out averages and stds 
        f_num_cells_spking_per_bin.write(str(avg_num_cells_spk[i]) + " ")
        f_num_cells_spking_per_bin.write(str(std_num_cells_spk[i]) + " ")
        f_num_cells_spking_per_intrabin.write(str(avg_num_cells_spk_intra[i]) + " ")
        f_num_cells_spking_per_intrabin.write(str(std_num_cells_spk_intra[i]) + " ")
        f_num_cells_spking_per_interbin.write(str(avg_num_cells_spk_inter[i]) + " ")
        f_num_cells_spking_per_interbin.write(str(std_num_cells_spk_inter[i]) + " ")
    
        f_num_spks_per_cell_per_bin.write(str(avg_num_spks_per_cell[i]) + " ")
        f_num_spks_per_cell_per_bin.write(str(std_num_spks_per_cell[i]) + " ")
        f_num_spks_per_cell_per_intrabin.write(str(avg_num_spks_per_cell_intra[i]) + " ")
        f_num_spks_per_cell_per_intrabin.write(str(std_num_spks_per_cell_intra[i]) + " ")
        f_num_spks_per_cell_per_interbin.write(str(avg_num_spks_per_cell_inter[i]) + " ")
        f_num_spks_per_cell_per_interbin.write(str(std_num_spks_per_cell_inter[i]) + " ")
    
        f_tot_num_spks_per_bin.write(str(avg_tot_num_spks[i]) + " ")
        f_tot_num_spks_per_bin.write(str(std_tot_num_spks[i]) + " ")
        f_tot_num_spks_per_intrabin.write(str(avg_tot_num_spks_intra[i]) + " ")
        f_tot_num_spks_per_intrabin.write(str(std_tot_num_spks_intra[i]) + " ")
        f_tot_num_spks_per_interbin.write(str(avg_tot_num_spks_inter[i]) + " ")
        f_tot_num_spks_per_interbin.write(str(std_tot_num_spks_inter[i]) + " ")
                  
        f_spk_hist_and_dist.write(' '.join(map(str,totspkhist_list[i])))
        f_spk_hist_and_dist.write("\n")
        f_spk_hist_and_dist.write(' '.join(map(str,totspkdist_smooth[i])))       
        
        f_num_cells_spking_per_bin.write(' '.join(map(str,num_cells_spk_per_bin[i]))) 
        f_num_cells_spking_per_intrabin.write(' '.join(map(str,num_cells_spk_per_intrabin[i])))
        f_num_cells_spking_per_interbin.write(' '.join(map(str,num_cells_spk_per_interbin[i])))
        
        f_num_spks_per_cell_per_bin.write(' '.join(map(str,num_spks_per_cell_per_bin[i]))) 
        f_num_spks_per_cell_per_intrabin.write(' '.join(map(str,num_spks_per_cell_per_intrabin[i])))
        f_num_spks_per_cell_per_interbin.write(' '.join(map(str,num_spks_per_cell_per_interbin[i])))
        
        f_tot_num_spks_per_bin.write(' '.join(map(str,tot_num_spks_per_bin[i])))
        f_tot_num_spks_per_intrabin.write(' '.join(map(str,tot_num_spks_per_intrabin[i])))
        f_tot_num_spks_per_interbin.write(' '.join(map(str,tot_num_spks_per_interbin[i])))
        
        f_FBF_per_bin.write(str(avg_fullburst_freq_per_bin_tot[i]) + " ")
        f_FBF_per_bin.write(str(std_fullburst_freq_per_bin_tot[i]) + " ")
        f_FBF_per_bin.write("\n")
        f_FBF_per_bin.write(' '.join(map(str,avg_FBF_per_bin_list[i])))
        f_FBF_per_bin.write("\n")
        f_FBF_per_bin.write(' '.join(map(str,std_FBF_per_bin_list[i])))
        
        f_intraBF_per_bin.write(str(avg_intraburst_freq_per_bin_tot[i]) + " ")
        f_intraBF_per_bin.write(str(std_intraburst_freq_per_bin_tot[i]) + " ")
        f_intraBF_per_bin.write("\n")
        f_intraBF_per_bin.write(' '.join(map(str,avg_intraBF_per_bin_list[i])))
        f_intraBF_per_bin.write("\n")
        f_intraBF_per_bin.write(' '.join(map(str,std_intraBF_per_bin_list[i])))
        
        f_interBF_per_bin.write(str(avg_interburst_freq_per_bin_tot[i]) + " ")
        f_interBF_per_bin.write(str(std_interburst_freq_per_bin_tot[i]) + " ")
        f_interBF_per_bin.write("\n")
        f_interBF_per_bin.write(' '.join(map(str,avg_interBF_per_bin_list[i])))
        f_interBF_per_bin.write("\n")
        f_interBF_per_bin.write(' '.join(map(str,std_interBF_per_bin_list[i])))
       
        f_num_cells_spking_per_bin.write("\n\n")
        f_num_spks_per_cell_per_bin.write("\n\n")
        f_tot_num_spks_per_bin.write("\n\n")
        f_FBF_per_bin.write("\n\n")
        f_spk_hist_and_dist.write("\n\n")
        
        f_num_cells_spking_per_intrabin.write("\n\n")
        f_num_spks_per_cell_per_intrabin.write("\n\n")
        f_tot_num_spks_per_intrabin.write("\n\n")
        f_intraBF_per_bin.write("\n\n")
               
        f_num_cells_spking_per_interbin.write("\n\n")
        f_num_spks_per_cell_per_interbin.write("\n\n")
        f_tot_num_spks_per_interbin.write("\n\n")
        f_interBF_per_bin.write("\n\n")
            
    f_num_cells_spking_per_bin.close()
    f_num_spks_per_cell_per_bin.close()
    f_tot_num_spks_per_bin.close()
    f_FBF_per_bin.close()
    f_spk_hist_and_dist.close()
    f_bins_ms.close()
    f_num_cells_spking_per_intrabin.close()
    f_num_spks_per_cell_per_intrabin.close()
    f_tot_num_spks_per_intrabin.close()
    f_intraBF_per_bin.close()
    f_intraburst_bins_ms.close()
    f_num_cells_spking_per_interbin.close()
    f_num_spks_per_cell_per_interbin.close()
    f_tot_num_spks_per_interbin.close()
    f_interBF_per_bin.close()
    f_interburst_bins_ms.close()
    
    ##write out values for each FBF/cell. 
    f_FBF_per_cell.write(' '.join(map(str,num_cells_w_FBF))) 
    f_FBF_per_cell.write("\n\n")                     
    for i in xrange(num_popln):
        f_FBF_per_cell.write(str(num_cells_w_FBF[i]) + " " + str(avg_fullburst_freq_per_cell_tot[i]) + " " + str(std_fullburst_freq_per_cell_tot[i]) + "\n")
    
    #f_FBF_per_cell.write("\n")    
    #for i in xrange(num_popln):
    #    for j in xrange(num_cells_w_FBF[i]):   #for each cell with at least one non-zero FBF
    #        f_FBF_per_cell.write(str(indices[i][j]) + " " + str((avg_FBF_per_cell_list[i][j])) + " " + str((std_FBF_per_cell_list[i][j])) + " ")   #write out average for jth cell
    #        for k in xrange(bin_ind_no_trans,numbins):    #bin number not including trans
    #            f_FBF_per_cell.write(str(fullburst_freq_nonzero[i][k,j])+" ")   #write out FBF for each bin (note:[k,j] since fullburst_freq has index as bin,cell  
    #        f_FBF_per_cell.write("\n")
    #    f_FBF_per_cell.write("\n\n")
    
    ##write out values for each intraBF/cell. 
    f_intraBF_per_cell.write(' '.join(map(str,num_cells_w_intraBF))) 
    f_intraBF_per_cell.write("\n\n")                     
    for i in xrange(num_popln):
        f_intraBF_per_cell.write(str(num_cells_w_intraBF[i]) + " " + str(avg_intraburst_freq_per_cell_tot[i]) + " " + str(std_intraburst_freq_per_cell_tot[i]) + "\n")
    
    #f_intraBF_per_cell.write("\n")    
    #for i in xrange(num_popln):
    #    for j in xrange(num_cells_w_intraBF[i]):   #for each cell with at least one non-zero FBF
    #        f_intraBF_per_cell.write(str(intraindices[i][j]) + " " + str((avg_intraBF_per_cell_list[i][j])) + " " + str((std_intraBF_per_cell_list[i][j])) + " ")   #write out average for jth cell
    #        for k in xrange(intrabin_ind_no_trans,num_intraburst_bins):    #bin number not including trans
    #            f_intraBF_per_cell.write(str(intraburst_freq_nonzero[i][k,j])+" ")   #write out FBF for each bin (note:[k,j] since fullburst_freq has index as bin,cell  
    #        f_intraBF_per_cell.write("\n")
    #    f_intraBF_per_cell.write("\n\n")
        
    ##write out values for each interBF/cell. 
    f_interBF_per_cell.write(' '.join(map(str,num_cells_w_interBF))) 
    f_interBF_per_cell.write("\n\n")                     
    for i in xrange(num_popln):
        f_interBF_per_cell.write(str(num_cells_w_interBF[i]) + " " + str(avg_interburst_freq_per_cell_tot[i]) + " " + str(std_interburst_freq_per_cell_tot[i]) + "\n")
    
    #f_interBF_per_cell.write("\n")    
    #for i in xrange(num_popln):
    #    for j in xrange(num_cells_w_interBF[i]):   #for each cell with at least one non-zero FBF
    #        f_interBF_per_cell.write(str(interindices[i][j]) + " " + str((avg_interBF_per_cell_list[i][j])) + " " + str((std_interBF_per_cell_list[i][j])) + " ")   #write out average for jth cell
    #        for k in xrange(interbin_ind_no_trans,num_interburst_bins):    #bin number not including trans
    #            f_interBF_per_cell.write(str(interburst_freq_nonzero[i][k,j])+" ")   #write out FBF for each bin (note:[k,j] since fullburst_freq has index as bin,cell  
    #        f_interBF_per_cell.write("\n")
    #    f_interBF_per_cell.write("\n\n")
        
    f_FBF_per_cell.close()
    f_intraBF_per_cell.close()
    f_interBF_per_cell.close()

