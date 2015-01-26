import csv
import os

configfile: "config.json"

COORDS = config["coord_file"]
REGION_NAME = config["region_name"]
NAME_MAPPING = config[config["dataset"]].get("name_mapping", "")

if config["data_type"] == "sunk":
    DATATYPE_TEXT = "sunk_"
elif config["data_type"] == "wssd":
    DATATYPE_TEXT = ""

REGRESS_NAME = "{dat}/{rn}_{type}{rn}.REGRESS.summary".format(dat = config["dataset"], rn = REGION_NAME, type = DATATYPE_TEXT)
GENOTYPE_FILE = config["genotype_file"]

TABLE_DIR = config["table_dir"]
PLOT_DIR = config["plot_dir"]

DIRS_TO_MAKE = ["log", TABLE_DIR, PLOT_DIR]

for folder in DIRS_TO_MAKE:
    if not os.path.exists(folder):
        os.makedirs(folder)

SCRIPT_DIR = config["script_dir"]
POP_FILE = config[config["dataset"]]["pop_file"]
GENE_GRAM_SETTINGS = config[config["dataset"]].get("gene_gram_settings", "")
SPP = config[config["dataset"]].get("spp", 500)

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
    input:  "%s/%s.table.tab" % (TABLE_DIR, REGION_NAME), 
            dynamic("%s/gene_grams/%s_%s_%s.{num_plot_type}" % (PLOT_DIR, REGION_NAME, config["dataset"], config["data_type"])),
            expand("%s/violins/%s_%s_violin_{name}_%s.{plot_type}" % (PLOT_DIR, REGION_NAME, config["dataset"], config["data_type"]),
            name = get_region_names(COORDS), plot_type = config["plot_type"]) 
    params: sge_opts=""

rule get_suns:
    input: GENOTYPE_FILE
    output: "%s/%s.table.tab" % (TABLE_DIR, REGION_NAME), temp("%s/%s.coords.tab" % (TABLE_DIR, REGION_NAME))
    params: sge_opts = "-l mfree=2G -N get_SUNs", suns = "/net/eichler/vol5/home/bnelsj/projects/gene_grams/hg19_suns.no_repeats_36bp_flanking.bed"
    shell:
        """sed -e '1d' {input} | cut -f 1-4 > {output[1]}
           bedtools intersect -a {output[1]} -b {params.suns} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
           awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' > {output[0]}
           sed -i '1ichr\tstart\tend\tname\tsize\tnSUNs' {output[0]}"""

rule plot_gene_grams:
    input: regions = GENOTYPE_FILE
    output: dynamic("%s/gene_grams/%s_%s_%s.{num_plot_type}" % (PLOT_DIR, REGION_NAME, config["dataset"], config["data_type"]))
    params: sge_opts = "-l mfree=8G -N gene_grams", output_prefix = "%s/gene_grams/%s_%s_%s" % (PLOT_DIR, REGION_NAME, config["dataset"], config["data_type"])
    run:
        for pt in config["plot_type"]:
            shell("""python {SCRIPT_DIR}/plotting/gene_gram.py {input.regions} {POP_FILE} {params.output_prefix} --plot_type {pt} --spp {SPP} {GENE_GRAM_SETTINGS}""")

rule plot_violins:
    input: "%s/%s_%s.genotypes.df" % (TABLE_DIR, REGION_NAME, config["data_type"])
    output: "%s/violins/%s_%s_violin_{name}_%s.{plot_type}" % (PLOT_DIR, REGION_NAME, config["dataset"], config["data_type"])
    params: sge_opts = "-l mfree=8G -N plot_violins"
    shell:
        """Rscript {SCRIPT_DIR}/plotting/genotype_violin.R {input} {output} {wildcards.name} {wildcards.plot_type} 3"""

rule get_long_table:
    input: regions = GENOTYPE_FILE
    output: "%s/%s_%s.genotypes.df" % (TABLE_DIR, REGION_NAME, config["data_type"])
    params: sge_opts = "-l mfree=8G -N make_long_table"
    shell:
        """Rscript {SCRIPT_DIR}/genotyper/transform_genotypes.R {input.regions} {POP_FILE} {output}"""

rule convert_genotypes:
    input: COORDS, REGRESS_NAME
    output: GENOTYPE_FILE, temp("%s_names.tab" % REGION_NAME), [temp(REGRESS_NAME + x) for x in [".tab", ".named.tab", ".header.tab"]]
    params: sge_opts="-l mfree=2G -N convert_genotypes"
    run:
        shell("""awk 'OFS="\t" {{ print $1":"$2"-"$3,$4 }}' {COORDS} | sort -k 1,1 > {REGION_NAME}_names.tab
                 sed 's/\s\+/\t/g' {REGRESS_NAME} | sed 1d | cut -f 5- | sort -k 1,1 > {REGRESS_NAME}.tab
                 join -j 1 {REGION_NAME}_names.tab {REGRESS_NAME}.tab | sed 's/\s\+/\t/g;s/[-:]/\t/g' > {REGRESS_NAME}.named.tab
                 head -n 1 {REGRESS_NAME} | sed 's/\s\+/\t/g' | cut -f 2-4,6- | awk 'OFS="\t" {{ $3=$3"\tname"; print }}' > {REGRESS_NAME}.header.tab
                 cat {REGRESS_NAME}.header.tab {REGRESS_NAME}.named.tab > {GENOTYPE_FILE}""")
        if NAME_MAPPING is not "":
            shell("""while read line; do set -- $line; sed -i "s/$1/$2/g" {GENOTYPE_FILE}; done < {NAME_MAPPING}""")
