FAMILY_NAME=$1

bedtools intersect -a $FAMILY_NAME.coords.bed -b /net/eichler/vol5/home/bnelsj/projects/gene_grams/hg19_suns.no_repeats_36bp_flanking.bed -wao | groupBy -g 1,2,3,4 -c 6,6 -o first,count | awk 'OFS="\t" {print $1, $2, $3, $4, $3-$2, $5 != "-1" ? $6 : 0}' > $FAMILY_NAME.table.tab; sed -i '1ichr\tstart\tend\tname\tsize\tnSUNs' $FAMILY_NAME.table.tab
