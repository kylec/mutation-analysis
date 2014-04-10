library(absCNseq)
library(DNAcopy)

# dna copy
args <- commandArgs(TRUE)

basename = args[1]
#basename = "FAP01-Vilar13-Vilar01"
varscanFile = paste(basename, "copynumber", sep=".")
dnacopyFile = paste(basename, "copynumber.dnacopy", sep=".")
cnFile = paste(basename, "abscn", sep=".")
snvFile = paste(basename, "mutect.keep.abscn_input", sep=".")
#purity,ploidy output
abscnFile = paste(basename, "absout", sep=".")

#if(FALSE) {
cat("read .called..")
cn <- read.table(varscanFile,header=T)
CNA.object <-CNA( genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file=dnacopyFile, row.names=F, col.names=F, quote=F, sep="\t")

cat("read dnacopy..")
# dnacopy to abscnseq output
d = segs2[,2:6]
colnames(d) = c("chrom", "start","stop", "num_seg", "seg_mean")

cat("sum bins..")
# sum bins
efflen =c()
for (i in 1:dim(d)[1]) {
  bins = cn[which(cn$chrom==d[i,]$chrom & cn$chr_start>=d[i,]$start & cn$chr_stop<=d[i,]$stop),]
  efflen = append(efflen, sum(bins$chr_stop-bins$chr_start))
}

#write output
ratio=2^d$seg_mean
#efflen = d$V3-d$V2
result = cbind(d[,1:3], efflen, ratio)
colnames(result) = c('chrom', 'loc.start', 'loc.end', 'eff.seg.len', 'normalized.ratio')
write.table(result, file=cnFile, row.names=F, quote=F, sep="\t")
#}

cat("run abscnseq..")
# run abscnseq
cat(abscnFile)
cat(snvFile)
my.res.list <- run.absCNSeq(cnFile, snvFile, basename, basename, seq.type="WES", min.seg.len=200)
write.table(my.res.list$searchRes, file=abscnFile, row.names=F, quote=F, sep="\t")

