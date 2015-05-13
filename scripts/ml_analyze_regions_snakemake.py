import sys
import socket
import os
import os.path
from optparse import OptionParser
from collections import defaultdict
import tempfile
import time
import subprocess
import gzip

import scipy as scp
import numpy as np
import tables
import pygr.Data

import scipy
import scipy.optimize
from scipy.special import gamma
from scipy.integrate import quadrature
from scipy.integrate import quad
from scipy.stats import norm,poisson

from math import *
import Bio.Statistics.lowess as biostats

from kitz_wssd.wssd_common_v2 import *

import ml_get_cpn as ml_methods
import genome_management.kg_file_handling as kgf

def file_exists(ls,file):
  for f in ls:
    if(f==file):
      return 1
  return 0

def mkdir(dir,file):

    ls_dir = os.listdir(dir)
    if(not(file_exists(ls_dir,file))):
        command = "mkdir %s/%s"%(dir,file)
        os.system(command)

    return "%s/%s"%(dir,file)

class ml_analyze_regions:
    
    def __init__(self,ml,in_locations,fn_genome_analysis_dir,fn_contigs,fn_mask,fn_out,combined_corrected_name,width,win_var=True,sunk_based=False,alt_wssd_contigs=None,alt_contig_map=None):
    
            fn_combined_adjusted_wssd = "%s/combined_corrected_wssd/%s"%(fn_genome_analysis_dir,combined_corrected_name)

            wssd_contigs = alt_wssd_contigs != None and alt_wssd_contigs or fn_contigs
            self.combined_adjusted_wssd = WssdFile(wssd_contigs,
                                                                                fnWssd=fn_combined_adjusted_wssd,
                                                                                overwrite=False,
                                                                                openMode='r')

            self.mask_wssd = DenseTrackSet(fn_contigs,
                                                                fnWssd=fn_mask,
                                                                overwrite=False,
                                                                openMode='r')

            #outdir = mkdir(fn_genome_analysis_dir,"ml_region_analysis")
            #FOUT = open("%s/%s"%(outdir,fn_out),'w')
            FOUT = open("%s"%(fn_out),'w')

            width = int(width)

            for location in in_locations:
                    if(location[0] == "#"):continue
                    #(TID,name,chr,start,end) = location.split()
                    #BED FILE
                    (chr,start,end) = location.split()[0:3]
                    

                    #print alt_contig_map
                    #print chr
                    #print self.mask_wssd.mContigNameLen
                    chr = alt_contig_map !=None and alt_contig_map[chr] or chr
                    TID = "%s:%s-%s"%(chr,start,end)
                    name = "%s:%s-%s"%(chr,start,end)
                    
                    if(not(chr in self.mask_wssd.mContigNameLen)):
                        print "skipping: ", name,chr,start,end
                        continue
                    print "loading region... subtracting 1 from start and end", name,chr,start,end
                    start = int(start)-1
                    end = min(int(end),int(self.mask_wssd.mContigNameLen[chr])-1)
                    #print "loading region... assuming half open 0 based coordinates", name,chr,start,end
                    #start = int(start)
                    #end = min(int(end),int(self.mask_wssd.mContigNameLen[chr]))
                    print start, end
                    
                    ######CHANGED IF SUNK    
                    if sunk_based:
                        corrected_depth = self.combined_adjusted_wssd.depth["wssd.combined_corrected"][chr][start:end,0,1]
                        masked_region = self.mask_wssd["isntSunkOrIsMasked"][chr][start:end]
                    else:
                        corrected_depth = self.combined_adjusted_wssd.depth["wssd.combined_corrected"][chr][start:end,:,0].sum(1)
                        masked_region = self.mask_wssd["mask"][chr][start:end,:].sum(1)>0
                    
                    #########BEFORE YOU FREAK OUT!!! READ THIS!!!! <RELAX> YOU'RE OK
                    #########YOU DID THE GENOTYPING CORRECTLY with adjustable windows... phew :)
                    if(win_var):
                        depth = corrected_depth[np.where(masked_region==0)]
                    else:
                        print "THIS IS WRONG!"
                        sys.exit(1)
                        depth = corrected_depth

                    print "finished loading"
                    wndstart = 0
                    len = depth.shape[0]
                    wndend = min(width,len)
                    
                    bywndcp =[]
                    regressed_bywndcp =[]
                    n_bases = []

                    n_windows = 0
                    w_ml_median = np.zeros(len)
                    w_regressed_median = np.zeros(len)
                    regressed_cp = ml.regressMu.regrinv(depth[wndstart:wndend].mean())

                    #################THIS IS THE FASTER CODE
                    if len!=0:    
                        depth_csum = np.cumsum(depth)         
                        idxs = np.r_[np.arange(min(width,len),len,width),len]-1
                        wnd_sizes = idxs-np.r_[0,idxs[0:idxs.shape[0]-1]]
                        c_wnd_depths = depth_csum[idxs]
                        wnd_depths = np.r_[c_wnd_depths[0],c_wnd_depths[1:]-c_wnd_depths[0:c_wnd_depths.shape[0]-1]]
                        wnd_means = wnd_depths/wnd_sizes
                        regwnds = ml.regressMu.regrinv(wnd_means)
                                
                        weighted_medians_idxs = np.cast['int32'](np.floor(np.arange(0,len)/float(width)))
                        #print len
                        #print weighted_medians_idxs
                        weighted_medians = regwnds[weighted_medians_idxs]
                            
                        #w_regressed_median = regwnds
                        #print "fast",np.median(w_regressed_median)
                        #print "fast",np.median(weighted_medians)
                        #print weighted_medians
                        w_regressed_median = weighted_medians
                        #print w_regressed_media$a
                    ##########################
                        
                    #########THIS IS THE ORIGINAL CODE
                    #if False:
                    #    while(wndstart<len-1):
                    #        wndend=min(len,wndstart+width)
                    #    
                        #    #if(wndend+width>len):
                    #    #    #    wndend = len
                    #    #    
                    #    #    #print wndstart,wndend
                    #    #    #cp,ll = ml.get_mle_cp(depth[wndstart:wndend])
                    #        cp,ll = 0,0
                    #        regressed_cp = ml.regressMu.regrinv(depth[wndstart:wndend].mean())
                    #    #    
                    #        regressed_bywndcp.append(float(regressed_cp))
                    #        bywndcp.append(float(cp))
                    #        n_bases.append(float(wndend-wndstart))
                    #    #
                    #        w_ml_median[wndstart:wndend] = float(cp)
                    #        w_regressed_median[wndstart:wndend] = float(regressed_cp)
                    #    #    
                    #        wndstart=wndend
                    #        n_windows+=1
                    ###################
                    
                    total_bases = np.array(n_bases).sum()
                    bywnd_cp = (np.array(bywndcp)*np.array(n_bases)/total_bases).sum()
                    bywnd_regressed_mean = (np.array(regressed_bywndcp)*np.array(n_bases)/total_bases).sum()
        
                    #print bywndcp
                    #print bywnd_regressed_mean
                    #print regressed_bywndcp
        
                    #if(n_windows==1):
                    #    bywnd_median = bywndcp[0]
                    #    bywnd_regressed_median = regressed_bywndcp[0]
                    #else:
                    #    bywnd_median = np.median(np.array(bywndcp))
                    #    bywnd_regressed_median = np.median(np.array(regressed_bywndcp))
                
                    #bywnd_median = np.median(w_ml_median) 
                    bywnd_median = 0 
                    bywnd_regressed_median = np.median(w_regressed_median)
                    #print w_regressed_median, name, chr 
                    #print bywnd_regressed_median
                    #print "3", w_regressed_median
            
                    #whole_reg_cp,ll = ml.get_mle_cp(depth)
                    whole_reg_cp,ll = 0,0
                    #whole_reg_regressed_cp = ml.regressMu.regrinv(depth.mean())
                    whole_reg_regressed_cp = 0
                    
                    #print bywndcp
                    #print bywnd_cp
                    #print round(bywnd_cp)
                    #print "%s %s %s %d %s %d %f %f %f %f %f %f\n"%(TID,name,chr,start+1,end,whole_reg_cp,bywnd_cp,bywnd_median,ll,whole_reg_regressed_cp,bywnd_regressed_mean,bywnd_regressed_median)
                    FOUT.write("\t".join(map(str,[TID,name,chr,start+1,end,whole_reg_cp,bywnd_cp,bywnd_median,ll,whole_reg_regressed_cp,bywnd_regressed_mean,bywnd_regressed_median])) + "\n")

            self.combined_adjusted_wssd.tbl.close()
            self.mask_wssd.tbl.close()


if __name__=='__main__':

    opts = OptionParser()
    opts.add_option('','--in_genomes',dest='fn_in_genomes')
    opts.add_option('','--in_regions',dest='fn_in_locations')
    opts.add_option('','--contigs',dest='fn_contigs')
    opts.add_option('','--mask_file',dest='fn_mask')
    opts.add_option('','--output_name',dest='fn_output_name')
    opts.add_option('','--sex_pop_index',dest='fn_sex_pop_index')
    opts.add_option('','--combined_corrected_name',dest='fn_comb_corr')
    opts.add_option('','--param_file',dest='fn_param')
    opts.add_option('','--output_dir',dest='out_dir',default=None)
    opts.add_option('','--omit_sample_name', dest='omit_sn', action='store_true', default=False)
    opts.add_option('','--window_width',dest='window_width',default=1000)
    opts.add_option('','--sunk_based',dest='sunk_based',action="store_true",default=False)
    opts.add_option('','--alt_primary_analysis_lambda',dest='alt_primary_analysis_lambda',default=None)
    
    opts.add_option('','--alt_wssd_contigs',dest='alt_wssd_contigs',default=None)
    opts.add_option('','--alt_contig_map',dest='alt_contig_map',default=None)
    #opts.add_option('','--finished',dest='fn_finished',default=None)
    
    (o, args) = opts.parse_args()

    genome_info = kgf.genome_info(o.fn_in_genomes,o.fn_sex_pop_index)

    in_genomes = open(o.fn_in_genomes,"r").readlines()
    in_locations = open(o.fn_in_locations,"r").readlines()
    alt_contig_map = o.alt_contig_map != None and dict([tuple([l.rstrip().split()[1],l.rstrip().split()[0]]) for l in open(o.alt_contig_map,'r').readlines()]) or None


    #primary_analysis_dir = "/net/gs/vol1/home/psudmant/ebod/psudmant/primary_analysis"
    
    #finished = []    
    #if(o.fn_finished !=None):
    #    finished_lines =open(o.fn_finished,'r').readlines()
    #    for line in finished_lines:
    #        genome_name = line.split("/")[9]
    #        print genome_name
    #        finished.append(genome_name)

    for in_genome in in_genomes:
        (genome,fn_wssd_dir,fn_bac_dir,chunk_dir,primary_analysis_dir) = in_genome.split()
        #if(o.fn_finished!=None):
        #    if(not(genome in finished)): 
        #        print "genome not complete yet",genome
        #        continue
        
        if o.alt_primary_analysis_lambda:
            fn_primary_analysis_mod = eval(o.alt_primary_analysis_lambda)
            primary_analysis_dir= fn_primary_analysis_mod(primary_analysis_dir)    


        if(genome_info.genomes[genome].passed_qc):
            genome_dir_name = genome.split(".")[0]

            wssd_lib_dir = "%s/%s"%(fn_wssd_dir,genome)
            bac_analysis_lib_dir = "%s/%s"%(fn_bac_dir,genome_dir_name)
            print genome, genome_dir_name
            #filtered_libs_dir = "%s/_filtered_libs_analysis.GC2"%(bac_analysis_lib_dir,fn_comb_corr)
            #if(not(os.path.exists(filtered_libs_dir))):
            #    print "error... filtered_libs_dir does not exist"
            fnfit_params = "%s/%s"%(bac_analysis_lib_dir,o.fn_param)
            #ml_class = ml_methods.ml_get_cp2_trunc_gaussian_twoPieceRegress(fnfit_params,lUseCp=[2,None],thruZero=True,max_cp=200)
        
            genome_analysis_dir = "%s/%s"%(primary_analysis_dir,genome)
            
            if o.out_dir == None:
                outdir = mkdir(fn_genome_analysis_dir,"ml_region_analysis")
                fn_out = "%s/%s"%(outdir,o.fn_output_name)
            else:
                if o.omit_sn:
                    fn_out = "%s/%s" % (o.out_dir,o.fn_output_name)
                else:
                    fn_out = "%s/%s_%s"%(o.out_dir,o.fn_output_name,genome)    
            
            ml_class = ml_methods.ml_get_cp2_trunc_gaussian(fnfit_params,lUseCp=[2,6],thruZero=True,max_cp=250)
            ml_analyze_regions(ml_class,in_locations,genome_analysis_dir,o.fn_contigs,o.fn_mask,fn_out,o.fn_comb_corr,width=o.window_width,sunk_based=o.sunk_based,alt_wssd_contigs=o.alt_wssd_contigs,alt_contig_map=alt_contig_map)
        else:
            print "FAILED"
