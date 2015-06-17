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
    parser.add_argument("--genotype_method", choices=["float", "GMM"], help="Output float or integer (Gaussian Mixture Model) genotypes")

    parser.add_argument("--subset", default=0)
    parser.add_argument("--total_subsets", default=1)

    parser.add_argument("--subset_indivs", nargs="+", help="Subset of individuals to genotype")
    parser.add_argumnet("--manifest", help="Path to manifest file with sample column")

    args = parser.parse_args()

#    (o, args) = opts.parse_args()
   
    max_cp = int(args.max_cp)
    subset = int(args.subset)
    total_subsets = int(args.total_subsets)

    tbx_dups = pysam.Tabixfile(args.fn_dup_tabix)
    GC_inf = GC_data(args.fn_GC_DTS, args.contig, args.fn_DTS_contigs)

    if args.subset_indivs is not None:
        indivs = args.subset_indivs
    elif args.manifest is not None:
        indivs = pd.read_table(args.manifest, header=0).sample.unique().tolist()
    else:
        indivs = list(pd.read_json("%s/gglob.idx" % args.gglob_dir).indivs)

    # GENOTYPE TIME!
    
    g = gt.genotyper(args.contig, gglob_dir = args.gglob_dir, plot_dir = args.plot_dir, subset_indivs = indivs, fn_fa=args.fn_fa, dup_tabix = tbx_dups, GC_inf = GC_inf)

    regions = pd.read_csv(args.fn_regions, header=None, delimiter="\t", index_col=None)
    regions.columns = ["chr", "start", "end", "name"]
    regions_by_contig = regions[regions['chr'] == args.contig]
    nregions = regions_by_contig.shape[0]

    FOUT = open(args.fn_out, 'w')
    if args.contig == args.header_chr and subset == 0:
        FOUT.write("chr\tstart\tend\tname\t%s\n"%("\t".join(indivs)))

    for i, row in regions_by_contig.iterrows():
        contig, s, e, name = row['chr'], int(row['start']), int(row['end']), row['name']
        if args.data_type == "wssd":
            X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        else:
            X, idx_s, idx_e = g.get_sunk_gt_matrix(contig, s, e)

        if args.genotype_method == "float":
            gt_list = np.mean(X, 1).tolist()
            gt_ordered = [gt_list[g.indivs.index(indiv)] for indiv in indivs]
            gts = "\t".join(map(str, gt_ordered))
        else:
            gts_by_indiv = g.simple_GMM_genotype(X, max_cp=max_cp)
            gts = "\t".join(["%d"%(gts_by_indiv[i]) for i in indivs])
        #gX.simple_plot("%s/%s_%d_%d.pdf"%(args.plot_dir, contig,s,e))
        FOUT.write("%s\t%d\t%d\t%s\t%s\n"%(contig, s, e, name, gts))
        print i, "%s\t%d\t%d\t%s\t%s\n"%(contig, s, e, name, gts)
