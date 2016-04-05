import csv
import os
import pandas as pd

configfile: "config.yaml"

shell.prefix("source dat/modules.txt; ")

COORDS = list(config.get("bedfiles").values())
REGION_NAMES = list(config.get("bedfiles").keys())

TABLE_DIR = config["table_dir"]
PLOT_DIR = config["plot_dir"]

REFERENCE = config["reference"]
DATASETS = config["datasets"]
DATATYPES = config["data_types"]

ds_manifest = pd.read_table(config["datasets_file"], header=0)
ds_manifest = ds_manifest.ix[ds_manifest.reference == REFERENCE, :]
ds_manifest.index = ds_manifest.dataset

DIRS_TO_MAKE = ["log", TABLE_DIR, PLOT_DIR]

genotyper = config.get("genotyper")

if genotyper == "wssd_cc":
    include: "workflows/wssd_cc_genotyper.snake"
elif genotyper == "gglob":
    include: "workflows/gglob_genotyper.snake"

for folder in DIRS_TO_MAKE:
    if not os.path.exists(folder):
        os.makedirs(folder)

GENE_GRAM_SETTINGS = config.get("gene_gram_settings", "")
SPP = config.get("spp", 500)

def get_region_names(region_names):
    names = []
    for region in region_names:
        with open(config.get("bedfiles")[region], 'r') as reader:
            for line in reader:
                names.append("%s.%s" % (region, line.rstrip().split()[3]))
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
                    return os.path.basename(coord_file).split(".")[0]

rule all:   
    input:  "%s/num_sunks.table.tab" % (TABLE_DIR),
            expand("%s/gene_grams/{fam}_{dataset}_{datatype}.0.{file_type}" % (PLOT_DIR),
            fam = REGION_NAMES, dataset = DATASETS, datatype = DATATYPES, file_type = config["plot_file_type"]),
            expand("%s/violin/{fam_name}.{dataset}_violin_{datatype}.{file_type}" % (PLOT_DIR),
            fam_name = get_region_names(REGION_NAMES), dataset = config["main_dataset"], datatype = DATATYPES, file_type = config["plot_file_type"]),
            expand("%s/{fam}.{plottype}_{datatype}.pdf" % PLOT_DIR, fam = REGION_NAMES, plottype=["violin", "scatter", "superpop"], datatype = DATATYPES),
            expand("{fam}/{fam}.{ds}.combined.{dt}.GMM.bed", fam = REGION_NAMES, ds = DATASETS, dt = DATATYPES)
    params: sge_opts=""

rule get_cn_wssd_variance:
    input:  gts = expand("{fam}/{dataset}/{dataset}.wssd.genotypes.tab", fam = REGION_NAMES, dataset = DATASETS),
            sunks = "%s/num_sunks.table.tab" % (TABLE_DIR)
    output: "%s/wssd_stats_by_family.tab" % (TABLE_DIR)
    params: sge_opts = "-l mfree=2G -N get_var", families = REGION_NAMES
    run:
        wssd_stats = pd.DataFrame(columns=["name"])
        datasets = DATASETS
        size = []
        for i, dataset in enumerate(datasets):
            wssd_mean, wssd_std, cn_two = [], [], []
            sample_dat = pd.read_csv(ds_manifest.loc[dataset]["manifest"], sep="\t", header=0)
            samples = sample_dat["sample"].tolist()
            for j, fam in enumerate(REGION_NAMES):
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
    input:  gts = expand("{fam}/{dataset}/{dataset}.sunk.genotypes.tab", fam = REGION_NAMES, dataset = DATASETS), 
            sunks = "%s/num_sunks.table.tab" % (TABLE_DIR)
    output: "%s/sunk_stats_by_region.tab" % TABLE_DIR
    params: sge_opts = "-l mfree=2G -N get_var", families = REGION_NAMES
    run:
        sunk_stats = pd.DataFrame(columns = ["name"])
        names = []
        for i, dataset in enumerate(DATASETS):
            sample_dat = pd.read_table(ds_manifest.loc[dataset]["manifest"], header=0)
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
        num_sunks = pd.read_csv(input.sunks, sep="\t", header=0)
        sunk_stats = num_sunks.merge(sunk_stats, on = "name")
        sunk_stats.to_csv(output[0], index=False, sep="\t")

rule get_sunks:
    input: COORDS
    output: "%s/num_sunks.table.tab" % (TABLE_DIR)
    params: sge_opts = "-l mfree=2G -N get_SUNKs", sunks = config["ref_files"][REFERENCE]["sunk_bed"]
    run:
        for i, coords in enumerate(COORDS):
            if i == 0:
                shell("""module load bedtools/2.21.0; bedtools intersect -a {coords} -b {params.sunks} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
                         awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' > {output[0]}""")
            else:
                shell("""module load bedtools/2.21.0; bedtools intersect -a {coords} -b {params.sunks} -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | 
                         awk 'OFS="\t" {{print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}}' >> {output[0]}""")
        shell("""sed -i '1ichr\tstart\tend\tname\tsize\tnSUNKs' {output[0]}""")

rule plot_gene_grams:
    input: expand("{fam}/{fam}.{dataset}.combined.{datatype}.bed", fam = REGION_NAMES, dataset = DATASETS, datatype = DATATYPES)
    output: "%s/gene_grams/{fam}_{dataset}_{datatype}.0.{file_type}" % (PLOT_DIR)
    params: sge_opts = "-l mfree=8G -N gene_grams"
    run:
        manifest = ds_manifest.loc[wildcards.dataset]["manifest"]
        shell("""python scripts/gene_gram.py {wildcards.fam}/{wildcards.fam}.{wildcards.dataset}.combined.{wildcards.datatype}.bed {manifest} {PLOT_DIR}/gene_grams/{wildcards.fam}_{wildcards.dataset}_{wildcards.datatype} --plot_type {wildcards.file_type} --spp {SPP} {GENE_GRAM_SETTINGS}""")

rule get_combined_pdfs:
    input: expand("%s/{fam}/violin_{datatype}.pdf" % PLOT_DIR, fam = REGION_NAMES, datatype = DATATYPES)
    params: ""

rule combine_violin_pdfs:
    input: expand("%s/{plottype}/{fam_name}.{dataset}_{plottype}_{datatype}.pdf" % (PLOT_DIR), plottype = ["violin", "scatter", "superpop"], fam_name = get_region_names(REGION_NAMES), dataset = config["main_dataset"], datatype = DATATYPES)
    output: "%s/{fam}.violin_{datatype}.pdf" % PLOT_DIR, "%s/{fam}.scatter_{datatype}.pdf" % PLOT_DIR, "%s/{fam}.superpop_{datatype}.pdf" % PLOT_DIR
    params: sge_opts = "-l mfree=8G -N pdf_combine"
    run:
        for pt in ["violin", "scatter", "superpop"]:
            outfile = [file for file in output if pt in file][0]
            shell("""gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={outfile} plots/{pt}/{wildcards.fam}*{pt}_{wildcards.datatype}.pdf""")

rule plot_violins:
    input: expand("%s/{fam}.{{dataset}}.{{datatype}}.genotypes.df" % (TABLE_DIR), fam = REGION_NAMES)
    output: "%s/violin/{fam_name}.{dataset}_violin_{datatype}.{file_type}" % (PLOT_DIR), "%s/scatter/{fam_name}.{dataset}_scatter_{datatype}.{file_type}" % (PLOT_DIR), "%s/superpop/{fam_name}.{dataset}_superpop_{datatype}.{file_type}" % (PLOT_DIR)
    params: sge_opts = "-l mfree=8G -N plot_violins"
    run:
        fam, name = wildcards.fam_name.split(".")[0], ".".join(wildcards.fam_name.split(".")[1:])
        input_table = [file for file in input if fam in file and wildcards.dataset in file][0]
        dat = pd.read_table(input_table)
        max_cp = max(dat.copy_num)
        max_cp = 7
        (coords, size) = get_coords_and_size_from_name(name, COORDS)
        title = "_".join([name, coords, size, config["reference"], wildcards.dataset, wildcards.datatype])
        shell("""Rscript scripts/genotype_violin.R {input_table} {output[0]} {name} {wildcards.file_type} {title} 3 violin super_pop_only --max_cp {max_cp} ; touch {output[0]}""")
        shell("""Rscript scripts/genotype_violin.R {input_table} {output[1]} {name} {wildcards.file_type} {title} 3 --max_cp {max_cp} ; touch {output[1]}""")
        shell("""Rscript scripts/genotype_violin.R {input_table} {output[2]} {name} {wildcards.file_type} {title} 3 super_pop_only --max_cp {max_cp} ; touch {output[2]}""")

rule get_tables:
    input: expand("%s/{fam}.{dataset}.{datatype}.genotypes.df" % (TABLE_DIR), fam = REGION_NAMES, dataset = DATASETS, datatype = DATATYPES)
    params: sge_opts = ""

rule get_long_table:
    input: regions = "{fam}/{fam}.{dataset}.combined.{datatype}.bed"
    output: "%s/{fam}.{dataset}.{datatype}.genotypes.df" % (TABLE_DIR)
    params: sge_opts = "-l mfree=8G -N make_long_table"
    run:
        master_manifest = config["master_manifest"]
        pop_codes = config["pop_codes"]
        shell("""Rscript scripts/transform_genotypes.R {input.regions} {master_manifest} {pop_codes} {wildcards.dataset} {output}""")

rule get_combined_GMM_genotypes:
    input: "{fam}/{fam}.{dataset}.combined.{datatype}.bed"
    output: "{fam}/{fam}.{dataset}.combined.{datatype}.GMM.bed"
    params: sge_opts = "-N GMM", max_cp = "12"
    shell:
        "python scripts/get_GMM_genotypes.py {input} {output} --max_cp {params.max_cp}"

rule combine_genotypes:
    input: expand("{{fam}}/{{fam}}.{ds}.{{datatype}}.genotypes.tab", ds = DATASETS)
    output: "{fam}/{fam}.{dataset}.combined.{datatype,\w+}.bed"
    params: sge_opts="-l mfree=2G -N combine_gt"
    run:
        fam = wildcards.fam
        dt = wildcards.datatype
        ds = wildcards.dataset
        fn_main = [x for x in input if ds in x][0]
        main_ds = pd.read_table(fn_main, na_values="NA")
        for app_ds in config["append_dataset"]:
            if app_ds != ds:
                append_dataset = pd.read_csv("{fam}/{fam}.{app_ds}.{dt}.genotypes.tab".format(fam=fam, app_ds=app_ds, dt=dt), header=0, sep="\t", index_col=False)
                main_ds = main_ds.merge(append_dataset, on=["chr", "start", "end", "name"])
        main_ds.to_csv("{fam}/{fam}.{ds}.combined.{dt}.bed".format(fam=fam, ds=ds, dt=dt), index=False, sep="\t", na_rep="NA")
