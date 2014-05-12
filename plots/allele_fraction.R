# plot allele fraction merge annotated

library(ggplot2)
library(reshape2)

# count loh after filtering
samples = system("ls *.ann.tsv | cut -d. -f1", inter=T)
for (i in 1:length(samples)) {
  basename = samples[i]
  a = read.table(paste(basename,"ann.tsv",sep="."), header=T, sep="\t")
  head(a)
  b = a[which(a$PN_B==1),]
  if (dim(b)[1] > 1) {
    mean_af = mean(b$AF_Tumor)
    sd_af = sd(b$AF_Tumor)
    #loh_af1 = mean_af - 2*sd_af
    #loh_af2 = mean_af + 2*sd_af
    sd1_count = dim(b[which(b$AF_Tumor <= mean_af-sd_af | b$AF_Tumor >= mean_af+sd_af),])[1]
    sd1_5_count = dim(b[which(b$AF_Tumor <= mean_af-1.5*sd_af | b$AF_Tumor >= mean_af+1.5*sd_af),])[1]
    sd2_count = dim(b[which(b$AF_Tumor <= mean_af-2*sd_af| b$AF_Tumor >= mean_af+2*sd_af),])[1]
    cat(paste(basename, round(mean_af,3), round(sd_af,3), dim(b)[1], sd1_count, sd1_5_count, sd2_count, "\n", sep="\t"))    
  }
}

# plot AF density/count
plot_list = c()
for (i in 1:length(samples)) {
  basename = samples[i]
  a = read.table(paste(basename,"ann.tsv",sep="."), header=T, sep="\t")
  head(a)
  b = a[which(a$PN_B==0 & a$f!="NA" & a$f <= 1),]
  b=b[grep("exonic", b$region_type),]
  p = ggplot(b) +
  geom_density(aes(x=AF_Tumor, fill="AF"), alpha=.3) + 
  geom_density(aes(x=f, fill="Adj. AF"), alpha=.3) + 
  #geom_histogram(aes(x=AF_Tumor), binwidth=.01, fill="blue", alpha=.2) +
  #geom_histogram(aes(x=f), binwidth=.01, fill="red", alpha=.2) +
  xlab("Allele_Fraction") + ggtitle(basename) + theme(text=element_text(size=14))
  plot_list[[i]] = p
  # using melt
  #c=cbind(b$AF_Tumor,b$f)
  #colnames(c)=c('a','b')
  #head(c)
  #d=melt(c)
  #ggplot(d, aes(x=value, fill=Var2)) + geom_density(alpha = .3)
}
library(grid)
library(gridExtra)
names(plot_list) = samples
png("test2.png", width=2000, height=2000)
arg_list <- c(plot_list, list(nrow=4, ncol=4))
do.call(grid.arrange, arg_list)
dev.off()