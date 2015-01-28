import csv
import os

configfile: "config.json"

COORDS = config["coord_file"]
REGION_NAME = config["region_name"]

TABLE_DIR = config["table_dir"]
PLOT_DIR = config["plot_dir"]

DIRS_TO_MAKE = [config["project_dir"], "log", TABLE_DIR, PLOT_DIR]

for folder in DIRS_TO_MAKE:
    if not os.path.exists(folder):
        os.makedirs(folder)

SCRIPT_DIR = config["script_dir"]
SPP = config.get("spp", 500)

def get_regress_name(wildcards, region_name = REGION_NAME):
    if wildcards.data_type == "sunk":
        datatype_text = "sunk_"
    elif wildcards.data_type == "wssd":
        datatype_text = ""
    return "{dat}/{rn}_{type}{rn}.REGRESS.summary".format(dat = wildcards.dataset, rn = region_name, type = datatype_text)

def get_region_names(coord_file):
    with open(coord_file, 'r') as reader:
        names = []
        for line in reader:
            names.append(line.rstrip().split()[3])
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
    input: expand("%s/%s_suns.table.tab" % (TABLE_DIR, REGION_NAME)), 
           expand("%s/violins/%s_{dataset}_violin_{name}_{data_type}.{plot_type}" % (PLOT_DIR, REGION_NAME), 
           dataset = config["datasets"], name = get_region_names(COORDS), data_type = config["data_types"], plot_type = config["plot_types"]),
           expand("%s/gene_grams/%s_{dataset}_{data_type}.0.{plot_type}" % (PLOT_DIR, REGION_NAME), 
           dataset = config["datasets"], data_type = config["data_types"], plot_type = config["plot_types"])
    params: sge_opts=""

rule get_suns:
    input: config["coord_file"] 
    output: "%s/%s_suns.table.tab" % (TABLE_DIR, REGION_NAME), temp("%s/%s.coords.tab" % (TABLE_DIR, REGION_NAME))
    params: sge_opts = "-l mfree=2G -N get_SUNs", suns = "/net/eichler/vol5/home/bnelsj/projects/gene_grams/hg19_suns.no_repeats_36bp_flanking.bed"
    shell:
        """cut -f 1-4 {input[0]} > {output[1]}
           bedtools intersect -a {output[1]} -b {params.suns} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
           awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' > {output[0]}
           sed -i '1ichr\tstart\tend\tname\tsize\tnSUNs' {output[0]}"""

rule plot_gene_grams:
    input: genotypes = "%s/%s_{dataset}_{data_type}_genotypes.tab" % (config["project_dir"], config["region_name"])
    output: "%s/gene_grams/%s_{dataset}_{data_type}.0.{plot_type}" % (PLOT_DIR, REGION_NAME)
    params: sge_opts = "-l mfree=8G -N gene_grams"
    run:
        pop_file = config[wildcards.dataset]["pop_file"]
        gene_gram_settings = config[wildcards.dataset]["gene_gram_settings"]
        output_prefix = "%s/gene_grams/%s_%s_%s" % (PLOT_DIR, REGION_NAME, wildcards.dataset, wildcards.data_type)
        shell("""python {SCRIPT_DIR}/plotting/gene_gram.py {input.genotypes} {pop_file} {output_prefix} --plot_type {wildcards.plot_type} --spp {SPP} {gene_gram_settings}""")

rule plot_violins:
    input: "%s/%s_{dataset}_{data_type}_genotypes.df" % (TABLE_DIR, REGION_NAME)
    output: "%s/violins/%s_{dataset}_violin_{name}_{data_type}.{plot_type}" % (PLOT_DIR, REGION_NAME)
    params: sge_opts = "-l mfree=8G -N plot_violins"
    shell:
        """Rscript {SCRIPT_DIR}/plotting/genotype_violin.R {input} {output} {wildcards.name} {wildcards.plot_type} 3"""

rule get_long_table:
    input: genotypes = "%s/%s_{dataset}_{data_type}_genotypes.tab" % (config["project_dir"], config["region_name"])
    output: "%s/%s_{dataset}_{data_type}_genotypes.df" % (TABLE_DIR, REGION_NAME)
    params: sge_opts = "-l mfree=8G -N make_long_table"
    run:
        pop_file = config[wildcards.dataset]["pop_file"]
        shell("""Rscript {SCRIPT_DIR}/genotyper/transform_genotypes.R {input.genotypes} {pop_file} {output}""")

rule convert_genotypes:
    input: coords = COORDS, regress = get_regress_name
    output: genotypes = "%s/%s_{dataset}_{data_type}_genotypes.tab" % (config["project_dir"], config["region_name"]), 
    params: sge_opts="-l mfree=2G -N convert_genotypes"
    run:
        name_mapping = config[wildcards.dataset].get("name_mapping", "")
        shell("""awk 'OFS="\t" {{ print $1":"$2"-"$3,$4 }}' {input.coords} | sort -k 1,1 > {REGION_NAME}_names.tab
                 sed 's/\s\+/\t/g' {input.regress} | sed 1d | cut -f 5- | sort -k 1,1 > {input.regress}.tab
                 join -j 1 {REGION_NAME}_names.tab {input.regress}.tab | sed 's/\s\+/\t/g;s/[-:]/\t/g' > {input.regress}.named.tab
                 head -n 1 {input.regress} | sed 's/\s\+/\t/g' | cut -f 2-4,6- | awk 'OFS="\t" {{ $3=$3"\tname"; print }}' > {input.regress}.header.tab
                 cat {input.regress}.header.tab {input.regress}.named.tab > {output.genotypes}""")
        if name_mapping is not "":
            shell("""while read line; do set -- $line; sed -i "s/$1/$2/g" {output.genotypes}; done < {name_mapping}""")
