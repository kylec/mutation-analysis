# plot mutation spectrum of fap

setwd("~/Projects/fap/mutect")

f = system("cat ~/Dropbox/lab_vilar/fap/samples_group.txt", intern=T)
patients_df = NULL
for (i in 1:length(f)) {
  samples = strsplit(f[i], "\t")[[1]]
  for (j in 1:length(samples)) {
    orig_sample = strsplit(samples[j], "_")[[1]][1]
    # shorten sample id
    sample = gsub("ilar", "", samples[j])
    fin = paste(orig_sample,"*keep",sep="")  
    mut_data = read.table(system(paste("ls ", fin), intern=T), header=T, sep="\t")
    
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
    
    sample_df = data.frame(type=type, fraction=fraction, sample=rep(sample, length(type)))
    patients_df = rbind(patients_df, sample_df)
  }
}

# set order of the type
patients_df$type = factor(patients_df$type, levels = c("C>T","C>A","C>G","T>C","T>A","T>G"))
names = levels(patients_df$sample)
ordered_names = c(names[grep("P", names)], names[grep("N", names)])
patients_df$sample = factor(patients_df$sample, levels=ordered_names)


# set color of bars
colors <- c("#F0E442", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#0072B2")
textsize = 10
p = ggplot(data=patients_df, aes(x=sample, y=fraction, fill=type, order=-as.numeric(type))) + geom_bar(stat="identity") +
  scale_fill_manual(values = colors) + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(breaks = seq(0,1,by=.1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size=textsize, angle=45, hjust=1),
        axis.text.y = element_text(size=textsize))

ggsave("fap_mut_spectrum.png", p)
