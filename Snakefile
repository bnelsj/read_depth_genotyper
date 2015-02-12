import csv
import os
import pandas as pd

configfile: "config.json"

COORDS = ["%s/%s.coords.bed" % (fam, fam) for fam in config["gene_families"]]

TABLE_DIR = config["table_dir"]
PLOT_DIR = config["plot_dir"]

DIRS_TO_MAKE = ["log", TABLE_DIR, PLOT_DIR]

for folder in DIRS_TO_MAKE:
    if not os.path.exists(folder):
        os.makedirs(folder)

SCRIPT_DIR = config["script_dir"]
GENE_GRAM_SETTINGS = config.get("gene_gram_settings", "")
SPP = config.get("spp", 500)

def get_pop_file(wildcards):
    return config[wildcards.dataset]["pop_file"]

def get_region_names(coord_files):
    names = []
    for coord_file in coord_files:
        with open(coord_file, 'r') as reader:
            for line in reader:
                names.append(line.rstrip().split()[3])
    return names

def get_coords_and_size_from_name(name, coord_files):
     for coord_file in coord_files:
        with open(coord_file, 'r') as reader:
            for line in reader:
                dat = line.rstrip().split()
                if dat[3] == name:
                    chr, start, end = dat[0], int(dat[1]), int(dat[2])
                    return ("%s:%d-%d" % (chr, start, end), str(end - start))

def get_family_from_name(name, coord_files):
    for coord_file in coord_files:
        with open(coord_file, 'r') as reader:
            for line in reader:
                test_name = line.rstrip().split()[3]
                if test_name == name:
                    return coord_file.replace(".coords.bed", "").split("/")[1]

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
    input:  "%s/num_suns.table.tab" % (TABLE_DIR),
            expand("%s/gene_grams/{fam}_{dataset}_{datatype}.0.{file_type}" % (PLOT_DIR),
            fam = config["gene_families"], dataset = config["dataset"], datatype = config["data_type"], file_type = config["plot_file_type"]),
            expand("%s/violin/{name}_{dataset}_violin_{datatype}.{file_type}" % (PLOT_DIR),
            dataset = config["dataset"], name = get_region_names(COORDS), datatype = config["data_type"], file_type = config["plot_file_type"]),
            expand("%s/{plottype}_{datatype}.pdf" % PLOT_DIR, plottype=["violin", "scatter"], datatype = config["data_type"])
    params: sge_opts=""

rule get_cn_wssd_variance:
    input:  gts = expand("{fam}/{dataset}/{dataset}_wssd_genotypes.tab", fam = config["gene_families"], dataset = config["dataset"]),
            suns = "%s/num_suns.table.tab" % (TABLE_DIR)
    output: "%s/wssd_stats_by_family.tab" % (TABLE_DIR)
    params: sge_opts = "-l mfree=2G -N get_var", families = config["gene_families"]
    run:
        wssd_stats = pd.DataFrame(columns=["name"])
        datasets = ["1kg", "hgdp"]
        size = []
        for i, dataset in enumerate(datasets):
            wssd_mean, wssd_std, cn_two = [], [], []
            sample_dat = pd.read_csv(config[dataset]["pop_file"], sep="\t", header=0)
            samples = sample_dat["sample"].tolist()
            for j, fam in enumerate(config["gene_families"]):
                data = pd.read_csv("%s/%s/%s_wssd_genotypes.tab" % (fam, dataset, dataset), sep="\t", header=0)
                common_samples = [sample for sample in samples if sample in data.columns]
                cns = data[common_samples]
                wssd_mean.append(cns.mean(axis=1).mean(axis=0))
                wssd_std.append(cns.std(axis=1).mean(axis=0))
                cn_two.append(cns.applymap(lambda x: x < 2.5).sum(axis=1).sum(axis=0))
                if i == 0:
                    size.append(int((data["end"] - data["start"]).mean()))
            wssd_stats[dataset + "_mean"] = wssd_mean
            wssd_stats[dataset + "_std"] = wssd_std
            wssd_stats[dataset + "_cn2"] = cn_two
        wssd_stats["name"] = list(params.families)
        wssd_stats["size"] = size
        wssd_stats.to_csv(output[0], index=False, sep="\t")

rule get_cn_sunk_variance:
    input:  gts = expand("{fam}/{dataset}/{dataset}_sunk_genotypes.tab", fam = config["gene_families"], dataset = config["dataset"]), 
            suns = "%s/num_suns.table.tab" % (TABLE_DIR)
    output: "%s/sunk_stats_by_region.tab" % TABLE_DIR
    params: sge_opts = "-l mfree=2G -N get_var", families = config["gene_families"]
    run:
        datasets = ["1kg", "hgdp"]
        sunk_stats = pd.DataFrame(columns = ["name"])
        names = []
        for i, dataset in enumerate(datasets):
            sample_dat = pd.read_csv(config[dataset]["pop_file"], sep="\t", header=0)
            samples = sample_dat["sample"].tolist()
            sunk_mean, sunk_std, cn_zero = [], [], []
            for j, fam in enumerate(params.families):
                data = pd.read_csv("%s/%s/%s_sunk_genotypes.tab" % (fam, dataset, dataset), sep="\t", header=0)
                common_samples = [sample for sample in samples if sample in data.columns]
                cns = data[common_samples]
                sunk_mean.extend(cns.mean(axis=1).tolist())
                sunk_std.extend(cns.std(axis=1).tolist())
                cn_zero.extend(cns.applymap(lambda x: x < 0.5).sum(axis=1).tolist())
                if i == 0:
                    names.extend(data["name"].tolist())
            sunk_stats[dataset + "_mean"] = sunk_mean
            sunk_stats[dataset + "_std"] = sunk_std
            sunk_stats[dataset + "_cn<0.5"] = cn_zero
        sunk_stats["name"] = names
        num_suns = pd.read_csv(input.suns, sep="\t", header=0)
        sunk_stats = num_suns.merge(sunk_stats, on = "name")
        sunk_stats.to_csv(output[0], index=False, sep="\t")

rule get_suns:
    input: expand("{fam}/{fam}.coords.bed", fam = config["gene_families"])
    output: "%s/num_suns.table.tab" % (TABLE_DIR)
    params: sge_opts = "-l mfree=2G -N get_SUNs", suns = "/net/eichler/vol5/home/bnelsj/projects/gene_grams/hg19_suns.no_repeats_36bp_flanking.bed"
    run:
        for i, fam in enumerate(config["gene_families"]):
            if i == 0:
                shell("""module load bedtools/2.21.0; bedtools intersect -a {fam}/{fam}.coords.bed -b {params.suns} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
                         awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' > {output[0]}""")
            else:
                shell("""module load bedtools/2.21.0; bedtools intersect -a {fam}/{fam}.coords.bed -b {params.suns} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
                         awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' >> {output[0]}""")
        shell("""sed -i '1ichr\tstart\tend\tname\tsize\tnSUNs' {output[0]}""")

rule plot_gene_grams:
    input: expand("{fam}/{fam}.{dataset}.combined.{datatype}.bed", fam = config["gene_families"], dataset = config["dataset"], datatype = config["data_type"])
    output: "%s/gene_grams/{fam}_{dataset}_{datatype}.0.{file_type}" % (PLOT_DIR)
    params: sge_opts = "-l mfree=8G -N gene_grams"
    run:
        pop_file = config[wildcards.dataset]["pop_file"]
        shell("""python {SCRIPT_DIR}/plotting/gene_gram.py {wildcards.fam}/{wildcards.fam}.{wildcards.dataset}.combined.{wildcards.datatype}.bed {pop_file} {PLOT_DIR}/gene_grams/{wildcards.fam}_{wildcards.dataset}_{wildcards.datatype} --plot_type {wildcards.file_type} --spp {SPP} {GENE_GRAM_SETTINGS}""")

rule get_combined_pdfs:
    input: expand("%s/violin_{datatype}.pdf" % PLOT_DIR, datatype = config["data_type"])
    params: ""

rule combine_violin_pdfs:
    input: expand("%s/{plottype}/{name}_{dataset}_{plottype}_{datatype}.pdf" % (PLOT_DIR), plottype = ["violin", "scatter"], name = get_region_names(COORDS), dataset = config["dataset"], datatype = config["data_type"])
    output: "%s/violin_{datatype}.pdf" % PLOT_DIR, "%s/scatter_{datatype}.pdf" % PLOT_DIR
    params: sge_opts = "-l mfree=8G -N pdf_combine"
    run:
        for pt in ["violin", "scatter"]:
            shell("""gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={PLOT_DIR}/{pt}_{wildcards.datatype}.pdf plots/{pt}/*{pt}_{wildcards.datatype}.pdf""")

rule plot_violins:
    input: expand("%s/{fam}_{{dataset}}_{{datatype}}.genotypes.df" % (TABLE_DIR), fam = config["gene_families"])
    output: "%s/violin/{name}_{dataset}_violin_{datatype}.{file_type}" % (PLOT_DIR), "%s/scatter/{name}_{dataset}_scatter_{datatype}.{file_type}" % (PLOT_DIR)
    params: sge_opts = "-l mfree=8G -N plot_violins"
    run:
        family = get_family_from_name(wildcards.name, COORDS)
        (coords, size) = get_coords_and_size_from_name(wildcards.name, COORDS)
        title = "_".join([wildcards.name, coords, size, config["build"], wildcards.dataset, wildcards.datatype])
        shell("""Rscript {SCRIPT_DIR}/plotting/genotype_violin.R {TABLE_DIR}/{family}_{wildcards.dataset}_{wildcards.datatype}.genotypes.df {output[0]} {wildcards.name} {wildcards.file_type} {title} 3 violin; touch {output[0]}""")
        shell("""Rscript {SCRIPT_DIR}/plotting/genotype_violin.R {TABLE_DIR}/{family}_{wildcards.dataset}_{wildcards.datatype}.genotypes.df {output[1]} {wildcards.name} {wildcards.file_type} {title} 3; touch {output[1]}""")

rule get_long_table:
    input: regions = "{fam}/{fam}.{dataset}.combined.{datatype}.bed"
    output: "%s/{fam}_{dataset}_{datatype}.genotypes.df" % (TABLE_DIR)
    params: sge_opts = "-l mfree=8G -N make_long_table"
    run:
        pop_file = config["master_manifest"]
        pop_codes = config["pop_codes"]
        shell("""Rscript {SCRIPT_DIR}/genotyper/transform_genotypes.R {input.regions} {pop_file} {pop_codes} {wildcards.dataset} {output}""")

rule convert_genotypes:
    input: "{fam}/{fam}.coords.bed"
    output: "{fam}/{fam}.{dataset}.combined.{datatype}.bed"
    params: sge_opts="-l mfree=2G -N convert_genotypes"
    run:
        fam = wildcards.fam
        datatype = wildcards.datatype
        for dataset in [wildcards.dataset, "archaics", "nhp"]:
            if datatype == "wssd":
                fn_gt = "_"
            else:
                fn_gt = "_sunk_"
            shell("""awk 'OFS="\t" {{ print $1":"$2"-"$3,$4 }}' {input[0]} | sort -k 1,1 > {fam}/{dataset}/{fam}_names.tab
                     sed 's/\s\+/\t/g' {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary | sed 1d | cut -f 5- | sort -k 1,1 > {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary.tab
                     join -j 1 {fam}/{dataset}/{fam}_names.tab {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary.tab | sed 's/\s\+/\t/g;s/[-:]/\t/g' > {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary.named.tab
                     head -n 1 {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary | sed 's/\s\+/\t/g' | cut -f 2-4,6- | awk 'OFS="\t" {{ $3=$3"\tname"; print }}' > {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary.header.tab
                     cat {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary.header.tab {fam}/{dataset}/{fam}{fn_gt}{fam}.REGRESS.summary.named.tab > {fam}/{dataset}/{dataset}_{datatype}_genotypes.tab""")
            name_mapping = config[dataset]["name_mapping"]
            if name_mapping != "":
                shell("""while read line; do set -- $line; sed -i "s/$1/$2/g" {fam}/{dataset}/{dataset}_{datatype}_genotypes.tab; done < {name_mapping}""")
            if dataset in ["1kg", "hgdp"]:
                shell("""cp {fam}/{dataset}/{dataset}_{datatype}_genotypes.tab {output[0]}.tmp""")
            else:
                shell("""cut -f 1-4 --complement {fam}/{dataset}/{dataset}_{datatype}_genotypes.tab | paste {output[0]}.tmp /dev/stdin > {output[0]}
                        cp {output[0]} {output[0]}.tmp""")
