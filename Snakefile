import csv
import os

configfile: "config.json"

if not os.path.exists(config["output_dir"]):
    os.makedirs(config["output_dir"])

SCRIPT_DIR = config["script_dir"]
POP_FILE = config[config["dataset"]]["pop_file"]
GENE_GRAM_SETTINGS = config[config["dataset"]]["gene_gram_settings"]
SPP = config[config["dataset"]]["spp"]

def get_region_names(region_file):
    with open(region_file, 'r') as reader:
        names = []
        csvreader = csv.DictReader(reader, delimiter="\t")
        for line in csvreader:
            names.append(line['name'])
    return names

def get_n_samples(region_file):
    with open(region_file, 'r') as reader:
        line = next(reader)
        nsamples = len(line.rstrip().split()) - 4
    return nsamples

def get_n_plots(region_file, spp):
    nsamples = get_n_samples(region_file)
    nplots = 0
    while nsamples > 0:
        nsamples -= spp
        nplots += 1
    return map(str, range(nplots))

rule all:
    input:  "%s/%s.table.tab" % (config["output_dir"], config["region_name"]), 
            expand("%s/%s_%s.{num}.{datatype}" % (config["output_dir"], config["region_name"], config["dataset"]), 
                    num = get_n_plots(config["region_file"], SPP), datatype = config["datatype"]),
            expand("%s/%s_%s_violin_{name}.{datatype}" % (config["output_dir"], config["region_name"], config["dataset"]), 
                    datatype = config["datatype"], name = get_region_names(config["region_file"]))
    params: sge_opts=""

rule get_suns:
    input: config["region_file"]
    output: "%s/%s.table.tab" % (config["output_dir"], config["region_name"]), temp("%s/%s.coords.tab" % (config["output_dir"], config["region_name"]))
    params: sge_opts = "-l mfree=2G -N get_SUNs", suns = "/net/eichler/vol5/home/bnelsj/projects/gene_grams/hg19_suns.no_repeats_36bp_flanking.bed"
    shell:
        """sed -e '1d' {input} | cut -f 1-4 > {output[1]}
           bedtools intersect -a {output[1]} -b {params.suns} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
           awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' > {output[0]}
           sed -i '1ichr\tstart\tend\tname\tsize\tnSUNs' {output[0]}"""

rule plot_gene_grams:
    input: regions = config["region_file"]
    output: "%s/%s_%s.{num}.{datatype}" % (config["output_dir"], config["region_name"], config["dataset"])
    params: sge_opts = "-l mfree=8G -N gene_grams", output_prefix = "%s/%s_%s" % (config["output_dir"], config["region_name"], config["dataset"])
    shell:
        """python {SCRIPT_DIR}/plotting/gene_gram.py {input.regions} {POP_FILE} {params.output_prefix} --plot_type {wildcards.datatype} {GENE_GRAM_SETTINGS}"""

rule plot_violins:
    input: "%s/%s.genotypes.df" % (config["output_dir"], config["region_name"])
    output: "%s/%s_%s_violin_{name}.{datatype}" % (config["output_dir"], config["region_name"], config["dataset"])
    params: sge_opts = "-l mfree=8G -N plot_violins"
    shell:
        """Rscript {SCRIPT_DIR}/plotting/genotype_violin.R {input} {output} {wildcards.name} {wildcards.datatype} 3"""

rule get_long_table:
    input: regions = config["region_file"]
    output: "%s/%s.genotypes.df" % (config["output_dir"], config["region_name"])
    params: sge_opts = "-l mfree=8G -N make_long_table"
    shell:
        """Rscript {SCRIPT_DIR}/genotyper/transform_genotypes.R {input.regions} {POP_FILE} {output}"""
