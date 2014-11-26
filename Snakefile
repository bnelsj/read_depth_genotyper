import os
import time

TOTAL_SUBSETS = 1
SCRIPT_DIR = "~bnelsj/src/read_depth_genotyper"
CHR_OUTPUT_DIR = "genotypes/chr"
FINAL_OUTPUT_DIR = "genotypes_all"
GGLOB_DIR = "/net/eichler/vol22/projects/1000_genomes_phase_II_III/nobackups/gglob"
REGIONS = "regions.bed"

MAX_CP = 50

CONTIGS = []
with open(REGIONS, 'r') as regions_file:
    for line in regions_file:
        dat = line.rstrip().split()
        chr = dat[0]
        if chr not in CONTIGS:
            CONTIGS.append(chr)

POP_FILE = "/net/eichler/vol2/eee_shared/1000_genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

DATA_TYPES = ["wssd", "sunk"]
GENOTYPE_METHODS = ["raw", "GMM"]

INIT_MODULES = "source %s/modules.txt" % SCRIPT_DIR

DIRS_TO_MAKE = ["log", CHR_OUTPUT_DIR, FINAL_OUTPUT_DIR]

for path in DIRS_TO_MAKE:
    if not os.path.exists(path):
        os.makedirs(path)

rule all:
    input: expand("%s/long/all.{method}.{type}.long" % FINAL_OUTPUT_DIR, method = GENOTYPE_METHODS, type = DATA_TYPES)
    params: sge_opts = ""

rule prep_for_plotting:
    input: expand("%s/long/all.{method}.{type}.long" % FINAL_OUTPUT_DIR, method = GENOTYPE_METHODS, type = DATA_TYPES)
    params: sge_opts = ""

rule make_long_tables:
    input: "%s/all.{method}.{type}" % FINAL_OUTPUT_DIR
    output: "%s/long/all.{method}.{type}.long" % FINAL_OUTPUT_DIR
    params: sge_opts = '-l mfree=8G -N long_tab'
    shell:
       "{INIT_MODULES}; Rscript {SCRIPT_DIR}/transform_genotypes.R {input} {POP_FILE} {output}"

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

rule genotype_by_chr:
    input: "%s/gglob.idx" % GGLOB_DIR
    output: "%s/{chr}.{method}.{type}" % CHR_OUTPUT_DIR
    params: sge_opts = "-l mfree=16G -N gt_{chr}_{method}_{type}", max_cp=str(MAX_CP), header_chr = CONTIGS[0]
    shell:
        "{INIT_MODULES}; python {SCRIPT_DIR}/combine_genotypes.py --regions {REGIONS} --contig {wildcards.chr} --output {output[0]} --gglob_dir {GGLOB_DIR} --genotype_method {wildcards.method} --data_type {wildcards.type} --max_cp {params.max_cp} --header_chr {params.header_chr}"
