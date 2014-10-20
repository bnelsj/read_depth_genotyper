import os

TOTAL_SUBSETS = 50
SCRIPT_DIR = "~bnelsj/genotyper_gaussian"
OUTPUT_DIR = "/net/eichler/vol22/projects/1000_genomes_phase_II_III/nobackups/bnelsj/1kg_genotypes"
GGLOB_DIR = "/net/eichler/vol22/projects/1000_genomes_phase_II_III/nobackups/jlhudd/dCGH/2014-02-05-1kg_vs_HGDP_redux/gglob"

CONTIGS = ['chr%s' % str(chr) for chr in list(range(1,23)) + ['X', 'Y']]

if not os.path.exists("log"):
    os.makedirs("log")

rule all:
    input: expand("%s/{chr}.genotypes" % OUTPUT_DIR, chr=CONTIGS)
    params: sge_opts = ""

rule combine_by_chr:
    input: expand("%s/{chr}.{s}.genotypes" % OUTPUT_DIR, chr=CONTIGS, s=list(range(TOTAL_SUBSETS)))
    output: "%s/{chr}.genotypes" % OUTPUT_DIR
    params: sge_opts = "-l mfree=4G"
    run:
        for fn_out in output:
            with open(fn_out, 'w') as out_genotype:
                files = []
                for infile in input:
                    if infile.split('/')[-1].split('.')[0] == "%s" % chr:
                        files.append(infile)
                    for file in files:
                        with open(file, 'r') as reader:
                            for line in reader:
                                out_genotype.write(line)

rule genotype_by_chr:
    input: 
    output: "%s/{chr}.{s}.genotypes" % OUTPUT_DIR
    params: sge_opts = "-l mfree=8G"
    shell:
        "python {SCRIPT_DIR}/combine_genotypes.py --contig {wildcards.chr} --output {output[0]} --gglob_dir {GGLOB_DIR} --subset {wildcards.s} --total_subsets {TOTAL_SUBSETS}"
