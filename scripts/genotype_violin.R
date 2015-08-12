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
region.name <- args[3]
output.type <- args[4]
plot.title <- args[5]

include_violin = "violin" %in% args
super_pop_only = "super_pop_only" %in% args

copy_nums <- read.table(input.file, header=TRUE)
filt.copy_nums <- copy_nums[copy_nums$name==region.name,]

legend_cn = filt.copy_nums
legend_cn$code = factor(legend_cn$code, levels = rev(unique(filt.copy_nums$code)), ordered = TRUE)

filt.copy_nums$code <- factor(filt.copy_nums$code, levels = rev(unique(filt.copy_nums$code)), ordered = TRUE)
sorted.copy_nums <- filt.copy_nums[with(filt.copy_nums, order(code, pop)),]
sorted.copy_nums$pop <- factor(sorted.copy_nums$pop, levels=unique(sorted.copy_nums$pop))
sorted.copy_nums$super_pop <- factor(sorted.copy_nums$super_pop, levels=unique(sorted.copy_nums$super_pop), ordered = TRUE)

code.order = rev(unique(filt.copy_nums$code))

sd.count <- sd(sorted.copy_nums$copy_num)
min.count <- floor(min(sorted.copy_nums$copy_num))
max.count <- ceiling(max(sorted.copy_nums$copy_num))
breaks = ifelse(max.count < 20, 0:max.count, seq(0, max.count, by=5))

threshold <- mean(filt.copy_nums$copy_num) + 8*sd(filt.copy_nums$copy_num)

xlab = "Super population"
ylab = "Copy number"

plot.legend <- ggplot(sorted.copy_nums, aes(x=code, y=copy_num, fill=code, name="Super population")) + 
        geom_violin() + theme_bw() + xlab(xlab) + ylab(ylab) + 
        scale_y_continuous(breaks=breaks, limits=c(-0.5, max.count+0.5)) + 
		guides(fill = guide_legend(reverse = TRUE))

pop.legend <- g_legend(plot.legend)

p2 <- ggplot(sorted.copy_nums, aes(x=pop, y=copy_num)) + 
    geom_point(alpha=0.5, colour='black', solid=T, size=1, position = position_jitter(h = 0, w=0.1)) + 
    theme_bw() + coord_flip() + xlab("Population") + ylab(ylab) + theme(legend.position="none", axis.text=element_text(size=6)) + 
    scale_y_continuous(breaks=breaks, limits=c(-0.5, max.count+0.5), minor_breaks=c())

if(include_violin) {
	p1 <- ggplot(sorted.copy_nums, aes(x=super_pop, y=copy_num, fill=code)) + 
    geom_violin() + geom_point(alpha=0.5, colour='black', solid=T, size=1, position = position_jitter(h = 0, w=0.1)) + 
    theme_bw() + coord_flip() + xlab(xlab) + ylab(ylab) + theme(legend.position="none") + 
    scale_y_continuous(breaks=breaks, limits=c(-0.5, max.count+0.5), minor_breaks=c())
} else {
	p1 <- ggplot(sorted.copy_nums, aes(x=super_pop, y=copy_num, fill=code)) + 
    geom_point(alpha=0.5, colour='black', solid=T, size=1, position = position_jitter(h = 0, w=0.1)) + 
    theme_bw() + coord_flip() + xlab(xlab) + ylab(ylab) + theme(legend.position="none") + 
    scale_y_continuous(breaks=breaks, limits=c(-0.5, max.count+0.5), minor_breaks=c())
}

if(super_pop_only) {
  if(output.type == "pdf") {
	  pdf(output.prefix, width=12, height=3)
	  grid.arrange(p1, pop.legend, ncol=2, widths=c(4/5, 1/5), main=plot.title)
	  dev.off()
  } else if(output.type == "png") {
  png(output.prefix, width=800, height=200)
  grid.arrange(p1, pop.legend, ncol=2, widths=c(4/5, 1/5), main=plot.title)
  dev.off()
  } else print(paste("Unsupported file type", output.type))
} else if(output.type == "pdf") {
    pdf(output.prefix, width=12, height=3*3)
    grid.arrange(arrangeGrob(p1, p2, ncol=2), pop.legend, ncol=2, widths=c(4/5, 1/5), main=plot.title)
    dev.off()
} else if(output.type == "png") {
png(output.prefix, width=800, height=200*3)
grid.arrange(arrangeGrob(p1, p2, ncol=2), pop.legend, ncol=2, widths=c(4/5, 1/5), main=plot.title)
dev.off()
} else print(paste("Unsupported file type", output.type))
