import argparse
import pandas as pd
import numpy as np
import simple_genotyper as gts


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genotypes", help="Tab-delimited input file with chr, start, end, name, and float cps")
    parser.add_argument("outfile")
    parser.add_argument("--max_cp", type = int, default = 12)

    args = parser.parse_args()

    dat = pd.read_table(args.genotypes)
    bed_header = ["chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"]
    of_bed_cols = ["chr", "start", "end", "name"]
    indivs = [col for col in dat.columns if col not in bed_header]

    with open(args.outfile, "w") as of:
        of.write("\t".join(of_bed_cols + indivs) + "\n")

        gs = gts.simple_genotyper()

        for i, row in dat.iterrows():
            X = map(float, row[indivs])
            GMM_cps = gs.simple_GMM_genotype(X, max_cp = args.max_cp)
            of_line = map(str, row[of_bed_cols].tolist() + GMM_cps)
            of.write("\t".join(of_line) + "\n")
