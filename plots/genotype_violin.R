library(ggplot2)
library(gridExtra)

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.prefix <- args[2]
region_name <- args[3]

copy_nums <- read.table(input.file, header=TRUE)
filt.copy_nums <- copy_nums[copy_nums$name==region_name,]
sorted.copy_nums <- filt.copy_nums[with(filt.copy_nums, order(super_pop, pop)),]
sorted.copy_nums$pop <- factor(sorted.copy_nums$pop, levels=unique(sorted.copy_nums$pop))

sd.count <- sd(sorted.copy_nums$copy_num)
min.count <- floor(min(sorted.copy_nums$copy_num))
max.count <- ceiling(max(sorted.copy_nums$copy_num))

threshold <- mean(filt.copy_nums$copy_num) + 8*sd(filt.copy_nums$copy_num)

xlab = "Super population"
ylab = "Copy number"

plot.legend <- ggplot(sorted.copy_nums, aes(x=super_pop, y=copy_num, fill=super_pop)) + geom_violin() + theme_bw() + coord_flip() + xlab(xlab) + ylab(ylab) + scale_y_continuous(breaks=0:max.count, limits=c(-0.5, max.count+0.5))

pop.legend <- g_legend(plot.legend)

p1 <- ggplot(sorted.copy_nums, aes(x=super_pop, y=copy_num, fill=super_pop)) + geom_violin() + geom_jitter(alpha=0.5, shape=1, position = position_jitter(h = 0, w=0.1)) + theme_bw() + coord_flip() + xlab(xlab) + ylab(ylab) + theme(legend.position="none") + scale_y_continuous(breaks=0:max.count, limits=c(-0.5, max.count+0.5))
p2 <- ggplot(sorted.copy_nums, aes(x=pop, y=copy_num, fill=super_pop)) + geom_violin() + theme_bw() + coord_flip() + xlab("Population") + ylab(ylab) + theme(legend.position="none") + scale_y_continuous(breaks=0:max.count, limits=c(-0.5, max.count+0.5))

pdf(paste0(output.prefix, ".pdf"), width=12, height=6)
grid.arrange(arrangeGrob(p1, p2, ncol=2), pop.legend, ncol=2, widths=c(5/6, 1/6), main=region_name)
dev.off()

png(paste0(output.prefix, ".png"), width=800, height=400)
grid.arrange(arrangeGrob(p1, p2, ncol=2), pop.legend, ncol=2, widths=c(5/6, 1/6), main=region_name)
dev.off()
