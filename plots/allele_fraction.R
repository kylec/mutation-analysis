# plot allele fraction merge annotated

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)


setwd("/Users/kchang3/Analysis/fap/expands/snv-only")
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

genes=c("KRAS", "APC", "GNAS", "TP53", "AKT1", "SOX9", "ARID1A", "CNOT3", "CDC27", "FBXW7", "TCF7L2", "GNAQ")
# plot AF density/count
plot_list = NULL
for (i in 1:length(samples)) {
  basename = samples[i]
  a = read.table(paste0(basename,".sps.ann.tsv"), header=T, sep="\t", stringsAsFactors=F)
  head(a)
  
  
  # allele frequency
  b = a
  b = b[grep("exonic", b$region_type),]
  # keep important gene names and plot them 
  b$region_name[b$region_name %in% genes == FALSE] = ""
  # plot gene names with AF_Tumor (allele fraction)
  p = ggplot(b) +
    geom_histogram(aes(x=AF_Tumor), binwidth=.01, fill="blue", alpha=.2) +
    geom_density(aes(x=AF_Tumor), alpha=.2) + 
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,15)) + 
    xlab("Allele_Fraction") + ggtitle(basename) + theme(text=element_text(size=14)) + 
    geom_text(aes(x=AF_Tumor, y=1, label=region_name, angle=90))
  p
  
  # absolute AF
  b = a
  b = b[which(b$PN_B==0 & b$f!="NA" & a$f <= 1),]
  b = b[grep("exonic", b$region_type),]
  b$region_name[b$region_name %in% genes == FALSE] = ""
  
  # number of clones
  sp = length(sort(unique(a$SP[!is.na(a$SP)])))
  
  p1 = ggplot(b) +
    geom_histogram(aes(x=f), binwidth=.01, fill="red", alpha=.2) +
    geom_density(aes(x=f), alpha=.2) + 
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,15)) + 
    xlab("Absolute Allele_Fraction") + ggtitle(basename) + theme(text=element_text(size=14)) + 
    geom_text(aes(x=f, y=1, label=region_name, angle=90))
  p1
   
  #g = arrangeGrob(p)
  g = arrangeGrob(p,p1)
  grid.arrange(g, ncol=1)
  # save image
  ggsave(paste0(basename, ".af.png"), g, width=10, height=15)
  
  #plot_list[[i]] = p
  # using melt
  #c=cbind(b$AF_Tumor,b$f)
  #colnames(c)=c('a','b')
  #head(c)
  #d=melt(c)
  #ggplot(d, aes(x=value, fill=Var2)) + geom_density(alpha = .3)
}

# plot all together
names(plot_list) = samples
png("test2.png", width=2000, height=2000)
arg_list <- c(plot_list, list(nrow=4, ncol=4))
do.call(grid.arrange, arg_list)
dev.off()
