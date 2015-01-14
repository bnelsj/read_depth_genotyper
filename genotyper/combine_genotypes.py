from optparse import OptionParser 
import argparse
import numpy as np
import pandas as pd
import csv
import pysam
import pdb
from bx.intervals.intersection import Interval, IntervalTree
import cluster
import genotyper as gt
from GC_data import GC_data

def get_regions_per_subset(regions, subset, total_subsets):
    subset_size = regions / total_subsets
    nremaining = regions - subset_size * total_subsets
    first_region = subset_size * (subset - 1) + min(nremaining, subset - 1)
    if subset <= nremaining:
        subset_size += 1
    last_region = first_region + subset_size - 1
    
    return (first_region, last_region)

def is_dup(row):
    gts = row.values[3:]
    gts = gts[gts!=-1]
    if np.sum(gts>2) != 0:
        return True
    return False


def bp(t, d = 1):
    return np.sum(t['end'].values - t['start'].values)/d

def overlap(s1, e1, s2, e2):
    l1 = e1 - s1
    l2 = e2 - s2
    mx_s = max(s1, s2)
    mn_e = min(e1, e2)
    o = float(mn_e - mx_s)
    of = min(o/l1, o/l2)
    return min(1.0, of)  


def merge_X_tables(X_M, X_F):
    
    X_files = [X_F, X_M]
    X_intervals = []
    X_tables = []
    for fn_X in X_files:
        i_tree = IntervalTree()
        t_gts = pd.read_csv(fn_X, header=0, delimiter="\t", index_col=None)
        X_tables.append(t_gts)
        for i, row in t_gts.iterrows():
            i_tree.insert_interval(Interval(row['start'], row['end'], i))
        X_intervals.append(i_tree) 
    
    calls = []
    
    added_from_table_2 = {}
    genotypes = []

    for i, row in X_tables[0].iterrows():
        s, e = row['start'], row['end']
        intersections = X_intervals[1].find(s,e)
        overlaps = np.array([overlap(s, e, interval.start, interval.end) for interval in intersections])
        #print s, e, overlaps, overlaps.shape[0]
        arg = overlaps.shape[0] == 0 and -1 or np.argmax(overlaps)
        if overlaps.shape[0]!=0 and overlaps[arg] >0.95:
            added_from_table_2[tuple([intersections[arg].start, intersections[arg].end])] = {'contig': "chrX", "start":s, "end":e}
        calls.append({'contig':"chrX", "start":s, "end":e})
            
    for i, row in X_tables[1].iterrows():
        s, e = row['start'], row['end']
        if not tuple([s,e]) in added_from_table_2:
            calls.append({'contig':"chrX", "start":s, "end":e})
    
    t = pd.DataFrame(calls)
    
    indivs = list(X_tables[0].columns[3:]) + list(X_tables[1].columns[3:])
    
    return t, indivs

def get_cps(t):
    indivs = t['indiv'].values
    inc = np.where(indivs != "WEA_Polish_ND15865_M") 
    cps =  t['cp'].values
    cps = cps[inc]
    indivs = indivs[inc]
    cps[np.isnan(cps)] = 0.0
    return cps, indivs

def genotype(gt, t, keystr, FOUT, FOUT_all, ordered_indivs):
    cps, indivs = get_cps(t) 
    cps = np.reshape(cps,(-1,1))
    gX = gt.GMM_genotype(cps, overload_indivs = indivs)
    #gX.simple_plot("./genotype_plots/%s.png"%keystr)
    
    mu = np.mean(cps)
    gts_by_indiv, gts_to_label, labels_to_gt = gX.get_gts_by_indiv()
    all_gts = np.array([gt for indiv, gt in gts_by_indiv.iteritems()])  
    
    n_0s = 0
    outstr = keystr

    if np.all(all_gts[0]==all_gts) and mu>=0.25: 
        for indiv in ordered_indivs:
            outstr = "%s\t2"%(outstr)
    else:
        for indiv in ordered_indivs:
            if gts_by_indiv[indiv] ==0: n_0s+=1
            outstr = "%s\t%d"%(outstr,gts_by_indiv[indiv])
    
    if n_0s!=len(indivs):
        FOUT.write("%s\n"%outstr)

    FOUT_all.write("%s\n"%outstr)
    #gt.output(FOUT, V_VCF, gX, 

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--contig", required=True)
    parser.add_argument("--output", dest="fn_out", required=True)
    parser.add_argument("--gglob_dir", required=True)
    parser.add_argument("--regions", dest="fn_regions", required=True)
    parser.add_argument("--plot_dir", default="plots")
    parser.add_argument("--fn_fa", default="/net/eichler/vol7/home/psudmant/genomes/fastas/hg19_1kg_phase2_reference/human_g1k_v37.fasta", help="reference genome fasta file (Default: %(default)s)")
    parser.add_argument("--GC_DTS", dest="fn_GC_DTS", default="/net/eichler/vol7/home/psudmant/genomes/GC_tracks/windowed_DTS/HG19/500_bp_slide_GC", help="GC tracks DTS file (Default: %(default)s")
    parser.add_argument("--DTS_contigs", dest='fn_DTS_contigs', default="/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/windowed_analysis/DTS_window_analysis/windows/hg19_slide/500_bp_windows.pkl.contigs", help="Contig sizes file (Default: %(default)s)")
    parser.add_argument("--dup_tabix", dest="fn_dup_tabix", default="/net/eichler/vol7/home/psudmant/genomes/annotations/hg19/superdups/superdups.merged.bed.gz", help="Superdups tabix file (Default: %(default)s)")
    parser.add_argument("--max_cp", default=12, type=int, help="Maximum cp to consider for GMM. Greater values will be rounded instead of fitted. Default: %(default)s")
    parser.add_argument("--header_chr", help="Name of chr to print header for")

    parser.add_argument("--data_type", choices=["wssd", "sunk"], help="Type of data to genotype (wssd or sunk)")
    parser.add_argument("--genotype_method", choices=["raw", "GMM"], help="Output raw (float) or integer (Gaussian Mixture Model) genotypes (choices: raw, GMM)")

    parser.add_argument("--subset", default=0)
    parser.add_argument("--total_subsets", default=1)

    parser.add_argument("--subset_indivs", nargs="+", help="Subset of individuals to genotype")

    args = parser.parse_args()

#    (o, args) = opts.parse_args()
   
    max_cp = int(args.max_cp)
    subset = int(args.subset)
    total_subsets = int(args.total_subsets)

    tbx_dups = pysam.Tabixfile(args.fn_dup_tabix)
    GC_inf = GC_data(args.fn_GC_DTS, args.contig, args.fn_DTS_contigs)

    if args.subset_indivs is None:
        indivs = list(pd.read_json("%s/gglob.idx" % args.gglob_dir).indivs)
    else:
        indivs = args.subset_indivs

    # GENOTYPE TIME!
    
    g = gt.genotyper(args.contig, gglob_dir = args.gglob_dir, plot_dir = args.plot_dir, subset_indivs = indivs, fn_fa=args.fn_fa, dup_tabix = tbx_dups, GC_inf = GC_inf)

    regions = pd.read_csv(args.fn_regions, header=None, delimiter="\t", index_col=None)
    regions.columns = ["chr", "start", "end", "name"]
    regions_by_contig = regions[regions['chr'] == args.contig]
    nregions = regions_by_contig.shape[0]

    FOUT = open(args.fn_out, 'w')
    if args.contig == args.header_chr and subset == 0:
        FOUT.write("contig\tstart\tend\tname\t%s\n"%("\t".join(indivs)))

    for i, row in regions_by_contig.iterrows():
        contig, s, e, name = row['contig'], int(row['start']), int(row['end']), row['name']
        if args.data_type == "wssd":
            X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        else:
            X, idx_s, idx_e = g.get_sunk_gt_matrix(contig, s, e)

        if args.genotype_method == "raw":
            gt_list = np.mean(X, 1).tolist()
            gt_ordered = [gt_list[g.indivs.index(indiv)] for indiv in indivs]
            gts = "\t".join(map(str, gt_ordered))
        else:
            gts_by_indiv = g.simple_GMM_genotype(X, max_cp=max_cp)
            gts = "\t".join(["%d"%(gts_by_indiv[i]) for i in indivs])
        #gX.simple_plot("%s/%s_%d_%d.pdf"%(args.plot_dir, contig,s,e))
        FOUT.write("%s\t%d\t%d\t%s\t%s\n"%(contig, s, e, name, gts))
        print i, "%s\t%d\t%d\t%s\t%s\n"%(contig, s, e, name, gts)
