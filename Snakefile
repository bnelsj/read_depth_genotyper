import os
import time

TOTAL_SUBSETS = 1
SCRIPT_DIR = "~bnelsj/genotyper_gaussian"
OUTPUT_DIR = "genotypes/subset"
CHR_OUTPUT_DIR = "genotypes/chr"
FINAL_OUTPUT_DIR = "genotypes_all"
GGLOB_DIR = "/net/eichler/vol22/projects/NCS_autism/nobackups/mduyzend/WGS_16p/gglob"
#REGIONS_FILE = "~bnelsj/gglob_get_genotypes/refGene.merged.bed"
REGIONS = "BOLA2.bed"
CONTIGS = ['chr%s' % str(chr) for chr in list(range(1,23)) + ['X', 'Y']]
HEADER_CHR = 'chr1'
def get_header(chr, HEADER_CHR, k):
    return "--header" if chr == HEADER_CHR and k == 0 else ''

DATA_TYPES = ["wssd", "sunk"]
GENOTYPE_METHODS = ["raw", "GMM"]

INIT_MODULES = "source ~bnelsj/genotyper_gaussian/modules.txt"

if not os.path.exists("log"):
    os.makedirs("log")

rule all:
    input: expand("%s/all.{method}.{type}" % FINAL_OUTPUT_DIR, method = GENOTYPE_METHODS, type = DATA_TYPES)
    params: sge_opts = ""

rule combine_by_chr:
    input: expand("%s/{chr}.{{method}}.{{type}}" % CHR_OUTPUT_DIR, chr=CONTIGS)
    output: "%s/all.{method}.{type}" % FINAL_OUTPUT_DIR
    params: sge_opts = '-l mfree=8G -N combine_all'
    run:
        with open(output[0], 'w') as outfile:
            for infile in input:
                with open(infile, 'r') as reader:
                    for line in reader:
                        outfile.write(line)
        time.sleep(10)

rule combine_by_subset:
    input: expand("%s/{{chr}}.{k}.{{method}}.{{type}}" % OUTPUT_DIR, k=list(range(TOTAL_SUBSETS)))
    output: "%s/{chr}.{method}.{type}" % CHR_OUTPUT_DIR
    params: sge_opts = "-l mfree=16G -N combine_{chr}"
    run:
        for fn_out in output:
            with open(fn_out, 'w') as out_genotype:
                files = []
                for infile in input:
                    with open(infile, 'r') as reader:
                        for line in reader:
                            out_genotype.write(line)
        time.sleep(10)

rule genotype_by_chr:
    input: "%s/gglob.idx" % GGLOB_DIR
    output: "%s/{chr}.{k}.{method}.{type}" % OUTPUT_DIR
    params: sge_opts = "-l mfree=16G -N gt_{chr}_{k}_{method}_{type}", max_cp="50", header=get_header("{chr}", HEADER_CHR, "{k}")
    shell:
        "{INIT_MODULES}; python {SCRIPT_DIR}/combine_genotypes.py --regions {REGIONS} --contig {wildcards.chr} --output {output[0]} --gglob_dir {GGLOB_DIR} --genotype_method {wildcards.method} --data_type {wildcards.type} --subset {wildcards.k} --total_subsets {TOTAL_SUBSETS} --max_cp {params.max_cp} {params.header}"
