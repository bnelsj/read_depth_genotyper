library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

genotypes.file <- args[1]
populations.file <- args[2]
pop.names.file <- args[3]
dataset <- args[4]
summary.file <- args[5]

# genotypes.file <- "all.raw.wssd"
# populations.file <- "/net/eichler/vol2/eee_shared/1000_genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
# summary.file <- "selected_loci_output.std_dataframe.txt"

data <- read.table(genotypes.file, header=TRUE, check.names=FALSE)
pop.names.table <- read.table(pop.names.file, header=TRUE, sep="\t")
populations <- read.table(populations.file, header=TRUE)

long.data <- melt(data, id.vars=c("chr", "start", "end", "name"), variable.name="sample", value.name="copy_num")
long.data <- long.data[with(long.data, order(sample)), ]

long.data$field = paste(long.data$chr, long.data$start, long.data$end, long.data$name, sep="_")

included.samples <- long.data$sample[long.data$sample %in% populations$sample]
long.data = long.data[long.data$sample %in% included.samples,]
#included.pops <- populations[populations$sample %in% included.samples,]
populations$code = as.character(populations$super_pop)

for(i in 1:length(row.names(populations))) {
  super_pop_i <- as.character(populations$code[i])
  if(dataset %in% pop.names.table$dataset) {
    codes = pop.names.table[pop.names.table$dataset == dataset,]
    if(super_pop_i %in% codes$super_pop) {
	  populations$code[i] = as.character(codes[codes$super_pop == super_pop_i, "code"])
    }
  }
}

populations$super_pop <- factor(populations$super_pop, levels = unique(populations$super_pop), ordered = TRUE)
populations <- populations[with(populations, order(populations$super_pop)), ]
#pops.merged <- merge(populations, included.pops, by="super_pop", all.x=TRUE, sort=FALSE)
pops.merged <- populations
merged.data <- merge(long.data, pops.merged, by="sample", all.x=TRUE, sort=FALSE)

std_dataframe <- merged.data[, c("super_pop", "code", "pop", "sample", "sex", "field", "name", "copy_num")]
std_dataframe$super_pop <- factor(std_dataframe$super_pop, levels = unique(pops.merged$super_pop), ordered = TRUE)
std_dataframe <- std_dataframe[with(std_dataframe, order(std_dataframe$super_pop)), ]

write.table(std_dataframe, file=summary.file, quote=FALSE, sep="\t", row.names=FALSE)
