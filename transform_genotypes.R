library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

genotypes.file <- args[1]
populations.file <- args[2]
summary.file <- args[3]

# genotypes.file <- "all.raw.wssd"
# populations.file <- "/net/eichler/vol2/eee_shared/1000_genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
# summary.file <- "selected_loci_output.std_dataframe.txt"

data <- read.table(genotypes.file, header=TRUE, check.names=FALSE)
long.data <- melt(data, id.vars=c("contig", "start", "end", "name"), variable.name="sample", value.name="copy_num")
long.data <- long.data[with(long.data, order(sample)), ]

long.data$field = paste(long.data$contig, long.data$start, long.data$end, long.data$name, sep="_")

length(long.data$sample)

populations <- read.table(populations.file, header=TRUE)
populations <- populations[with(populations, order(sample)), ]
length(populations$sample)

merged.data <- merge(long.data, populations, by="sample", all.x=TRUE)
length(merged.data$sample)

std_dataframe <- merged.data[, c("super_pop", "pop", "sample", "gender", "field", "copy_num")]
length(std_dataframe$sample)

write.table(std_dataframe, file=summary.file, quote=FALSE, sep="\t", row.names=FALSE)
