INPUT_FILE=$1
OUTPUT_PREFIX=$2

bedtools intersect -a $INPUT_FILE -b /net/eichler/vol5/home/bnelsj/projects/gene_grams/hg19_suns.no_repeats_36bp_flanking.bed -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | awk 'OFS="\t" {print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}' > $OUTPUT_PREFIX.tmp.tab
bedtools intersect -a $INPUT_FILE -b /net/eichler/vol2/eee_shared/assemblies/hg19/sunks/hg19_sunks.bed -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | awk 'OFS="\t" {print $5 != "-1" ? $6 : 0}' | paste $OUTPUT_PREFIX.tmp.tab - > $OUTPUT_PREFIX.table.tab
sed -i '1ichr\tstart\tend\tname\tsize\tnSUNs\tnSUNKs' $OUTPUT_PREFIX.table.tab
rm $OUTPUT_PREFIX.tmp.tab
