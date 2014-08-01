# Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'

# Get a list of the bedtools output files you'd like to read in

## @knitr readplot
dir="/Users/kchang3/Analysis/fap_ampliseq/"
setwd(paste0(dir,"cov"))
files <- list.files(pattern="all.txt$")

# Optional, create short sample names from the filenames. 
# For example, in this experiment, my sample filenames might look like this:
# prefixToTrash-01.pe.on.pos.dedup.realigned.recalibrated.bam
# prefixToTrash-02.pe.on.pos.dedup.realigned.recalibrated.bam
# prefixToTrash-03.pe.on.pos.dedup.realigned.recalibrated.bam
# This regular expression leaves me with "samp01", "samp02", and "samp03" in the legend.
labs <- paste("samp", gsub("prefixToTrash-0|\\.pe\\.on\\.pos\\.dedup\\.realigned\\.recalibrated\\.bam\\.cov\\.hist\\.txt\\.all\\.txt", "", files, perl=TRUE), sep="")

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i], sep="\t", quote="")
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}


# Pick some colors
# Ugly:
# cols <- 1:length(cov)
# Prettier:
# ?colorRampPalette
# display.brewer.all()
#library(RColorBrewer)
#cols <- brewer.pal(length(cov), "Dark2")
cols = rainbow(52)
# Save the graph to a file
#png("exome-coverage-plots.png", h=1000, w=1000, pointsize=20)

# Create plot area, but do not plot anything. Add gridlines and axis labels.
thr=4000
thr1=thr+1

plot(cov[[1]][2:thr1, 2], cov_cumul[[1]][1:thr], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage")
abline(v = 100, col = "gray60")
abline(v = 200, col = "gray60")
abline(v = 500, col = "gray60")
abline(v = 800, col = "gray60")
abline(v = 1000, col = "gray60")
abline(v = 2000, col = "gray60")
abline(v = thr, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
#s(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:thr1, 2], cov_cumul[[i]][1:thr], type='l', lwd=3, col=cols[i])

# Add a legend using the nice sample labeles rather than the full filenames.
#legend("topright", legend=labs, col=cols, lty=1, lwd=4)
#dev.off()

## @knitr summary
sample_info=read.table(paste0(dir,"samples.txt"), header=T, sep="\t", quote="")
samples=gsub('.bam.hist.all.txt', '', files)
dd=NULL
for (i in 1:length(cov)) {
    cov_sum = summary(rep(cov[[i]]$V2, cov[[i]]$V3))
    type=sample_info[sample_info$sample_id==samples[i], ]$type
    #print(c(samples[i], as.character(type), cov_sum[3], cov_sum[4]), sep="\t", quote=FALSE)
    dd = rbind(dd, c(samples[i], as.character(type), as.numeric(cov_sum[3]), as.numeric(cov_sum[4])))
}
ddd = data.frame(dd)
colnames(ddd) = c("sample", "type", "median", "mean")
ddd$median = as.numeric(as.character(ddd$median))
ddd$mean = as.numeric(as.character(ddd$mean))
library(ggplot2)
textsize=9
dddord=ddd[order(ddd$median), ]
dddord$sample = factor(dddord$sample, levels = dddord$sample)
#ggplot(data=dddord, aes(x=sample, y=median, fill=type)) + geom_point(aes(colour=type)) +  
#  theme(legend.position="top", axis.text.x = element_text(size=textsize, angle=45, hjust=1),
#        axis.text.y = element_text(size=textsize))

#ggplot(data=dddord, aes(x=sample, y=mean, fill=type)) + geom_point(aes(colour=type)) +  
#    theme(legend.position="top", axis.text.x = element_text(size=textsize, angle=45, hjust=1),
#          axis.text.y = element_text(size=textsize))

library(reshape2)
mdd = melt(dddord, id=c("sample","type"))
ggplot(data=mdd, aes(x=sample, y=value, fill=variable)) + geom_point(aes(colour=variable)) +  
    geom_abline(intercept=c(500,1000)) + 
    theme(legend.position="top", axis.text.x = element_text(size=textsize, angle=45, hjust=1),
          axis.text.y = element_text(size=textsize))


## @knitr mutect
setwd(paste0(dir,"iot-caller"))
a= read.table("cov.txt", header=T)
head(a)
hist(a$af, breaks=seq(0,1,by=.005), xaxt='n',main="Allele fraction of somatic snv", xlab="Alelle fraction")
axis(side=1, at=seq(0,1,by=.05))

