# plot gene heatmap
# Kyle Chang

library(ggplot2)
library(reshape2)
arg = commandArgs()
file = arg[0]
#file ="/Users/kyle_air/Dropbox/hgsc/gene.txt"
data = read.table(file, sep="\t", header=T)

# ggplot2 plots based on the order of factor levels in your variables. When a file
# is loaded into data.frame, factor level of "sample" is sorted by default, we can
# change the order of sample by setting it as the original order or other order desired

# default sorted factor level of sample when file is loaded
# head(data$sample)
#[1] AA-3555 AA-3979 AA-3986 AA-3681 AG-3611 AA-A00W
#91 Levels: A6-2672 A6-2674 A6-2676 A6-2678 A6-2683 A6-3808 AA-3516 AA-3518 AA-3519 ... AG-A01L

# change factor level to my intended order (mutually exclusive)
data$sample = factor(data$sample, levels = data$sample)
#head(data$sample)
#[1] AA-3555 AA-3979 AA-3986 AA-3681 AG-3611 AA-A00W
#91 Levels: AA-3555 AA-3979 AA-3986 AA-3681 AG-3611 AA-A00W AG-3598 AA-3852 AA-3522 ... AA-3562

# format data into sample, gene, and mutation_status columns
format_data = melt(data)

png("ordered_gene_heatmap.png", width=1200, height=400)
# text size
base_size <- 12
(p <- ggplot(format_data, aes(sample, variable)) + geom_tile(aes(fill = value),colour = "black") 
 + scale_fill_gradient(low = "white",high = "steelblue") 
 + theme_grey(base_size = base_size) 
 + labs(x = "", y= "") + scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) 
 + theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.7, angle = 90, hjust = 0, colour = "grey50")))  

dev.off()
