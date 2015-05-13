import sys
import socket
import os
import os.path
from optparse import OptionParser
#import scipy as scp
import numpy as np
import matplotlib.pyplot as plt
import pylab

import genome_management.kg_file_handling as kgf
import math

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

class region_info:
    def __init__(self,name,chr,start,end,TID):
        self.name = name
        self.chr = chr
        self.start = start
        self.end = end
        self.frequencies_by_pop = {}
        self.cps_by_genome = {}
        self.transcript_id = TID
        self.TID = TID
        self.cps_all = []
        self.pop_by_genome = {}

    def add_info_from_genome(self,cp,genome):
        if(not(genome.pop in self.frequencies_by_pop)):
            self.frequencies_by_pop[genome.pop] = []
        self.frequencies_by_pop[genome.pop].append(cp)
        self.cps_by_genome[genome.genome_name] = cp
        self.pop_by_genome[genome.genome_name] = genome.pop
        self.cps_all.append(cp)


  #def get_var(self):
   # self.vars = {}
        #self.cps_all = np.array(self.cps_all) 
  #  varT = self.cps_all.var()
    
#        self.vars["all"]=varT
  
#    self.means = {}                 
 #   meanT = self.cps_all.mean(1)
  #  self.means["all"] = meanT       
                                    
   # for pop,copies_by_pop in self.frequencies_by_pop.iteritems():
        #    copies_by_pop = np.array(copies_by_pop)
     # self.vars[pop] = self.summary[:,pop_index].var(1)
     # self.means[pop] = self.summary[:,pop_index].mean(1)
    

#    self.vsts = {}
 #   self.fsts = {}
  #  for pop,pop_index in self.indivs_by_pop.iteritems():
   #   for pop_2,pop_index_2 in self.indivs_by_pop.iteritems():
    #    n_pop = float(pop_index.shape[0])
     #   n_pop_2 = float(pop_index_2.shape[0])
      #  both_pops = np.r_[self.indivs_by_pop[pop],self.indivs_by_pop[pop_2]]
       # var_both = self.summary[:,both_pops].var(1)
       # N = n_pop+n_pop_2
       # self.vsts["_".join([pop,pop_2])] = (var_both - ((self.vars[pop]*n_pop+self.vars[pop_2]*n_pop_2)/N)) / var_both


def make_output_file(region,region_info,outdir,cell_line_info,genome_info):
    outfile_name = "%s/%s_pop_summary.csv"%(outdir,region_info.name)

    FOUT = open(outfile_name,'w')
    FOUT.write("indiv,cp,pop,cell lines fixed, cell lines in Nitrogen,coverage\n")

    for indiv,cp in region_info.cps_by_genome.iteritems():
        pop = region_info.pop_by_genome[indiv]
        output = indiv in cell_line_info and cell_line_info[indiv] or ""
        output = "%s,%d,%s,%s,%f\n"%(indiv,cp,pop,output,genome_info.genomes[indiv].coverage)
        FOUT.write(output)
        print output    

def make_simple_plot(region,region_info,outdir,cell_line_info,genome_info):

    plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
    f=plt.figure()
    f.set_figwidth(6)
    f.set_figheight(6)

    axescolor  = '#f6f6f6'
    
    left, width = 0.1, 0.8
    rect1 = [left, 0.1, width, 0.8] #left, bottom, width, height
    ax = f.add_axes(rect1)

    colors = {'Yoruba':'r','European':'b','Asian':'g'}

    for indiv,cp in region_info.cps_by_genome.iteritems():
        cvg = genome_info.genomes[indiv].coverage
        fixed_cell_line = cell_line_info[indiv].split(",")[0].rstrip() == "yes"
        liquid_nitrogen_cell_line = cell_line_info[indiv].split(",")[1].rstrip() == "yes"
        color = colors[genome_info.genomes[indiv].pop]
        ax.plot(np.array([cvg]),np.array([cp]),'%so'%(color))

    ax.set_xlabel("cvg",size=20)
    ax.set_ylabel("copy",size=20)
    ax.set_title("%s"%(region_info.name),size=20)
    f.savefig("%s/%s_copy_vs_cvg.pdf"%(outdir,region_info.name),format='pdf')    
    plt.close(1)

def make_histogram(region,region_info,outdir,great_ape_gene_hashes):

    print region_info.name

    plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
    f=plt.figure()
    f.set_figwidth(10)
    f.set_figheight(10)
    nbins=0
    mx=0
    mn=100
    
    do_apes=True

    great_ape_cps = {}
    if do_apes:
        for ape,gene_hash in great_ape_gene_hashes.iteritems():
            if not region_info.TID in gene_hash:
                do_apes=False
                print "ID does not exist for APE"
                print region_info.TID
                break
            great_ape_cps[ape] = gene_hash[region_info.TID]
            mx=int(max(great_ape_cps[ape],mx))
            mn=int(min(great_ape_cps[ape],mn))
    
    axescolor  = '#f6f6f6'
    
    left, width = 0.1, 0.8
    rect1 = [left, 0.1, width, 0.8] #left, bottom, width, height

    for pop,freq_info in region_info.frequencies_by_pop.iteritems():
        #nbins = int(round(max(nbins,max(freq_info))))
        mx=int(max(max(freq_info),mx))
        mn=int(min(min(freq_info),mn))

    #nbins+=1
    nbins = mx-mn+1

    labels = []

    pop_to_hists = {}
    for pop,freq_info in region_info.frequencies_by_pop.iteritems():
        print pop,freq_info
        pop_to_hists[pop] = np.histogram(np.array(freq_info),bins=nbins,range=[mn,mx],normed=True,new=True)[0]
        print  np.histogram(np.array(freq_info),bins=nbins,range=[mn,mx],normed=True,new=True)
        print pop_to_hists[pop]

    x = np.arange(mn,mx+1)
    width=.25
    print x
    
    for i in range(x.shape[0]):
        labels.append(str(x[i]))    
    
    ax = f.add_axes(rect1)
    bars = {}
    
    leg = []
    leg_colors = []
    lines = []

    k=0
    colors = ['r','g','b','o']
    
    starty = .9
    sub=.03
    i=0
    for pop,freqs in region_info.frequencies_by_pop.iteritems():
        med = np.median(np.array(freqs))
        sig2 = np.array(freqs).var()
        leg.append("%s med: %d var: %.1f"%(pop,int(med),sig2))
        i+=1

    for pop,hist in pop_to_hists.iteritems():
        bars[pop] = ax.bar(x+k*width,hist,width,color=colors[k],alpha=0.5)
        leg_colors.append(colors[k])
        #ax.legend(bars[pop][0],pop)
        lines.append(bars[pop][0])
        k+=1

    ape_colors = ['orange','purple','yellow','brown']
    k=0
    if do_apes:
        for ape,cp in great_ape_cps.iteritems():
            bars_ape = ax.bar(np.array([cp]),np.array([.1]),width/2,color=ape_colors[k],alpha=.8)
            leg.append("%s %f"%(ape,cp))
            lines.append(bars_ape[0])
            k+=1

    ax.set_xticks(x+width*k/2)
    ax.set_xticklabels(labels,size=20)

    ax.grid(color='k',linestyle='--',linewidth=1,alpha=.3)    
    yticklabels = [str(x) for x in np.arange(0,1,.1)]
    ax.set_yticklabels(yticklabels,size=20)
    ax.set_ylabel("%",size=20)
    ax.set_xlabel("cp number",size=20)
    ax.legend(lines,leg)
    ax.set_title("%s"%(region_info.name),size=20)

    f.savefig("%s/%s_pop_hist.pdf"%(outdir,region_info.name),format='pdf')    
    plt.close(1)

    return

    k=0
    
    for pop,ihist in percent_hists.iteritems():
        percent_hists[pop] = ihist/ihist.sum()
        #jhplot(x,hist,"|%s"%(colors[k]))
        #hist(x)
        vlines(x+float(k)/3,zeros,percent_hists[pop],color=colors[k],linewidth=7)
        k+=1
        leg.append(pop)
    
    #legend(leg)
    title("percent")
    print leg
    legend(leg)
    
    f.get_axes()[0].xaxis.set_ticks(range(21))
    #f.add_axes([0,40,0,1],xticks=[0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20],label='axis2',axisbg='g')
    #[0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20])
    
    f=figure(2)
    k=0

    for pop,ihist in mode_hists.iteritems():
        mode_hists[pop] = ihist/ihist.sum()
        #plot(x,hist,"|%s"%(colors[k]))
        #hist(x)
        vlines(x+float(k)/5,zeros,mode_hists[pop],color=colors[k],linewidth=7)
        k+=1
    legend(leg)
    title("Predicted copy number %s"%(name))
    xlabel("predicted copy number")
    ylabel("percentage of population")
    f.get_axes()[0].xaxis.set_ticks(range(21))
    savefig("%smode_hist.png"%(name),format='png')    


    print percent_hists
    print mode_hists

def load_plot_regions(fn_regions):
    if fn_regions == None: return []
    plot_regions = []
    for line in open(fn_regions,'r').readlines():
        if line[0] == "#": continue
        print line
        sline = line.split()
        uID = "%s:"%(sline[1])
        uID += ":".join(sline[2:5])
        plot_regions.append(uID)
        print uID

    return plot_regions

def get_transcript_ids(fn_transcript_id):

    print fn_transcript_id
    gene_id_list = open(fn_transcript_id,'r').readlines()
    transcript_ids = {} 
    for gene_info in gene_id_list:
        (TID,name,chr,start,end,unmasked_len,GCp) = gene_info.split()        
        transcript_ids["%s:%s:%s"%(chr,start,end)] = {"tid":TID,"chr":chr,"start":start,"end":end,"unmasked":unmasked_len,"GC":GCp}

    return transcript_ids


def get_cp_by_gene(gene_file):
    
    cps_by_TID = {}
    for line in open(gene_file,'r').readlines():
        if len(line.split()) == 0: continue
        (chr,start,end,TID,cp) = line.split()            
        cps_by_TID[TID] = float(cp)

    return cps_by_TID

def get_calkan_cp_calls(fn_great_ape_cps_files):

    calkan_cp_calls = {}

    if(fn_great_ape_cps_files!=None):
        for line in open(fn_great_ape_cps_files,'r').readlines():
            (genome,gene_file) = line.split()                    
            calkan_cp_calls[genome] = get_cp_by_gene(gene_file)

    return calkan_cp_calls


if __name__=='__main__':

    opts = OptionParser()
    opts.add_option('','--input_file_name',dest='input_file_name')
    opts.add_option('','--input_genomes',dest='fn_input_genomes')
    opts.add_option('','--outdir',dest='outdir')
    opts.add_option('','--sex_pop_index',dest='fn_sex_pop_index')
    #opts.add_option('','--analysis_dir',dest='fn_analysis_dir')
    opts.add_option('','--input_regions',dest='input_regions',default=None)
    opts.add_option('','--out_file',dest='outfile',default=None)
    opts.add_option('','--regress',dest='regress',action='store_true',default=False)
    opts.add_option('','--plot_regions',dest='plot_regions',default=None)
    opts.add_option('','--do_plotting',action="store_true",dest='do_plotting',default=False)
    opts.add_option('','--great_ape_cps_files',dest='fn_great_ape_cps_files',default=None)
    opts.add_option('','--cell_line_information',dest='fn_cell_line_info',default=None)
    opts.add_option('','--output_coverage',dest='output_cvg',action='store_true',default=False)
    opts.add_option('','--simple_plot',dest='simple_plot',action='store_true',default=False)
    opts.add_option('','--input_dir',dest='input_dir',default=None)


    #opts.add_option('','--transcript_id_file',dest='fn_transcript_id')
    #opts.add_option('','--call_metric',dest='outfile',default="summary")

    #opts.add_option('','--out_genomes',dest='fn_out_genomes')
    (o, args) = opts.parse_args()

    great_ape_cps = get_calkan_cp_calls(o.fn_great_ape_cps_files)


    cell_line_info = {}
    if o.fn_cell_line_info != None:
        read_cell_line_info = open(o.fn_cell_line_info,'r').readlines()
        for cell_line_line in read_cell_line_info:
            (name,cells_fixed,in_nitrogen) = cell_line_line.split(",")                            
            cell_line_info[name] = "%s,%s"%(cells_fixed,in_nitrogen.rstrip())
            print cell_line_info[name]

    mkdir("./",o.outdir)
    print "loading genome information"
    genome_info = kgf.genome_info(o.fn_input_genomes,o.fn_sex_pop_index,QC_check=o.output_cvg)
    print "done"

    regions_by_uID = {}
    #print o.input_regions

    expected_len = 0
    if o.input_regions != None:
        for l in open(o.input_regions,'r').readlines():
            expected_len+= (l[0]!="#") and 1 

    input_genomes = open(o.fn_input_genomes,'r').readlines()
    plot_regions = load_plot_regions(o.plot_regions)

    outstr = "\t".join(["name", "chr", "start", "end", "TID"])
    for input_genomes_line in input_genomes:
        (genome_id,fn_wssd_dir,fn_bac_dir,chunk_dir,primary_analysis_dir) = input_genomes_line.split()    
        if genome_id[0] == "#": continue
        genome_ob = genome_info.genomes[genome_id]

        if o.input_dir is None:
            input_file = "%s/%s/ml_region_analysis/%s"%(primary_analysis_dir,genome_id,o.input_file_name)
        else:
            input_file = "%s/%s_%s"%(o.input_dir,o.input_file_name,genome_id)

        print input_file

        ##########check the output file exists
        #if(not(os.path.exists("%s/%s/ml_region_analysis/%s"%(primary_analysis_dir,genome_id,o.input_file_name)))):
        if not os.path.exists(input_file):
            print "%s does not appear to exist" % (input_file)
            print     
            print '%s my have failed previous QC or may still be running' % (genome_id)
            continue
        
        ##############check the output file is of the correct length
        #################here we coudl also put "take the first n"
        #analyzed_by_ml_lines = open("%s/%s/ml_region_analysis/%s"%(primary_analysis_dir,genome_id,o.input_file_name)).readlines()
        analyzed_by_ml_lines = open(input_file, "r").readlines()
        if(len(analyzed_by_ml_lines) != expected_len):
            print "expected:%d encountered:%d" % (expected_len, len(analyzed_by_ml_lines))
            print "expected number of lines in %s does not match that in %s" % (analyzed_by_ml_lines, o.input_regions)
            #continue
        
        print "\t getting information %s" %(genome_id)
        outstr += "\t%s" % genome_id
        
        for analysis_line in analyzed_by_ml_lines:
            (name,TID,chr,start,end,cp,bywnd_cp,median,ll,regressed_cp,regressed_cp_by_wnd,regressed_cp_median) = analysis_line.split()
            if o.regress:
                cp = float(regressed_cp_median)
            else:
                cp = float(median)
            
            uID = "%s:%s:%s:%s"%(TID,chr,start,end)
            if(not(uID in regions_by_uID)):
                regions_by_uID[uID] = region_info(name,chr,start,end,TID)
            regions_by_uID[uID].add_info_from_genome(cp,genome_ob)                

    outstr+="\n"
    for region_uID, region_inf in regions_by_uID.iteritems():
        
        outstr+="\t".join([region_inf.name,region_inf.chr,region_inf.start,region_inf.end,region_inf.transcript_id])
        #for genome_id,genome in genome_info.genomes.iteritems():
        for input_genomes_line in input_genomes:
            (genome_id,fn_wssd_dir,fn_bac_dir,chunk_dir,primary_analysis_dir) = input_genomes_line.split()    
            if genome_id[0] =="#": continue
            if genome_id in region_inf.cps_by_genome:
                #print genome_id
                outstr+="\t%f"%(region_inf.cps_by_genome[genome_id])
            else:
                print "ERROR genome_id not in region_info"
                print genome_id
                print region_inf.cps_by_genome
                sys.exit(1)
        outstr+="\n"    

#    print outstr
    if o.outfile != None: 
        open("%s/%s"%(o.outdir,o.outfile),'w').write(outstr)
        #print percent_hists[pop]
        #print hist
#        percent_hists[pop]=ihist + percent_hists[pop]
#        mode_hists[pop][np.where(ihist==np.amax(ihist))[0]]+=1 


