# plot mutation spectrum of fap samples
library
setwd("~/Projects/fap/mutect")

# evs coordinates to keep.
keep_coords = read.table("/Users/kyle_air/Projects/fap/expands/fap_mutect_novel_01_Jun03_155401.report", header=T, sep="\t", fill=T)
keep_coords = keep_coords[,c("chr", "hg19_pos")]
# mutect files extension
input_ext = "*.keep.exome" 

f = system("cat ~/Dropbox/lab_vilar/fap/samples_group.txt", intern=T)
patients_df = NULL
novel_dbsnp_count_df = NULL

for (i in 1:length(f)) {
  samples = strsplit(f[i], "\t")[[1]]
  for (j in 1:length(samples)) {
    orig_sample = strsplit(samples[j], "_")[[1]][1]
    # shorten sample id
    sample = gsub("ilar", "", samples[j])
    fin = paste(orig_sample, input_ext,sep="")  
    mut_data = read.table(system(paste("ls ", fin), intern=T), header=T, sep="\t")
    mut_data[,"contig"] = gsub("chr", "", mut_data[,"contig"])
    
    #filter mutations by evs
    mut_data = merge(mut_data, keep_coords, by.x=c("contig", "position"), by.y=c("chr", "hg19_pos"))
    
    mut_spect = paste(mut_data$ref_allele, mut_data$alt_allele, sep=">")
    count(mut_spect)
    mut_spect = gsub("A>C", "T>G", mut_spect)
    mut_spect = gsub("A>G", "T>C", mut_spect)
    mut_spect = gsub("A>T", "T>A", mut_spect)
    mut_spect = gsub("G>A", "C>T", mut_spect)
    mut_spect = gsub("G>C", "C>G", mut_spect)
    mut_spect = gsub("G>T", "C>A", mut_spect)
    
    # data frame of mutspect, counts, sample names
    count = count(mut_spect)$freq
    type = count(mut_spect)$x
    fraction = count/sum(count)
    patients_df = rbind(patients_df, data.frame(type=type, fraction=fraction, count=count, sample=rep(sample, length(type))))
    
    # cosmic/dbsnp/novel count
    novel = count(mut_data$dbsnp_site)
    novel_dbsnp_count_df = rbind(novel_dbsnp_count_df, data.frame(sample=rep(sample, dim(novel)[1]), type=novel$x, count=novel$freq))
    
  }
}

# set order of the type
patients_df$type = factor(patients_df$type, levels = c("C>T","C>A","C>G","T>C","T>A","T>G"))
names = levels(patients_df$sample)
ordered_names = c(names[grep("P", names)], names[grep("N", names)])
patients_df$sample = factor(patients_df$sample, levels=ordered_names)
novel_dbsnp_count_df$sample = factor(novel_dbsnp_count_df$sample, levels=ordered_names)

# set color of bars
colors <- c("#F0E442", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#0072B2")
textsize = 10
mutspect_plot = ggplot(data=patients_df, aes(x=sample, y=fraction, fill=type, order=-as.numeric(type))) + geom_bar(stat="identity") +
  scale_fill_manual(values = colors) +  
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(breaks = seq(0,1,by=.1), expand = c(0, 0)) +
  theme(legend.position="top", axis.text.x = element_text(size=textsize, angle=45, hjust=1),
        axis.text.y = element_text(size=textsize))

mutcount_plot = ggplot(data=novel_dbsnp_count_df, aes(x=sample, y=count, fill=type)) + geom_bar(stat="identity") +
  scale_fill_brewer(palette="Set1") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(breaks=seq(0,y,by=25), expand = c(0, 0)) +
  theme(legend.position="top", axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=textsize))

g = arrangeGrob(mutcount_plot, mutspect_plot)
grid.arrange(g, ncol=1)

ggsave("fap_mutation_spectrum.png", g)
