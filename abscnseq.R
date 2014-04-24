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

# read raw copynumber
cn <- read.table(varscanFile,header=T)

# segment copynumber
cat("segment copynumber..")
if(file.exists(dnacopyFile)) {
  segs = read.table(dnacopyFile, sep="\t")
  colnames(segs) = c("chrom", "start","stop", "num_seg", "seg_mean")
} else {
  CNA.object <-CNA( genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
  CNA.smoothed <- smooth.CNA(CNA.object)
  segs <- segment(CNA.smoothed, verbose=0, min.width=2)
  segs2 = segs$output
  write.table(segs2[,2:6], file=dnacopyFile, row.names=F, col.names=F, quote=F, sep="\t")
  d = segs2[,2:6]
  colnames(d) = c("chrom", "start","stop", "num_seg", "seg_mean")
}

# lei's code
# dnacopy to abscnseq cnv input
cat("convert segmented data..")

rd.list <- split(cn,cn[,"chrom"])
segs.list <- split(segs,segs[,"chrom"])
cat("loop...")
tmp.list <- list()
n.windows <- 0
for(i in 1:length(rd.list)) {
  rd.df <- rd.list[[i]]
  seg.df <- segs.list[[i]]
  n.seg <- nrow(seg.df)
  
  seg.len <- integer(n.seg)
  seg.r <- numeric(n.seg)
  for(s in 1:n.seg) {
    
    seg.rd <- rd.df[rd.df[,2]>=seg.df[s,2] & rd.df[,2]<=seg.df[s,3],]  # all start location
    
    n.windows <- n.windows+nrow(seg.rd)
    
    seg.len[s] <- sum(seg.rd[,4]) # effective length without counting gaps
    wts <- seg.rd[,4]/seg.len[s]
    tmp.r <- seg.rd[,6]/seg.rd[,5]
    seg.r[s] <- sum(wts*tmp.r) # normalized ratio
    
  }
  
  tmp.list[[i]] <- cbind(seg.df[,1:3],seg.len, seg.r)
}

seg.data <- do.call(rbind,tmp.list)
colnames(seg.data) <- c("chrom", "loc.start", "loc.end", "eff.seg.len", "normalized.ratio")
write.table(seg.data, file=cnFile, row.names=F, quote=F, sep="\t")

cat("run abscnseq..")
# run abscnseq
cat(abscnFile)
cat(snvFile)
my.res.list <- run.absCNSeq(cnFile, snvFile, basename, basename, seq.type="WES", min.seg.len=200)
write.table(my.res.list$searchRes, file=abscnFile, row.names=F, quote=F, sep="\t")
#plot.absCN(seg.CN, chromnum=1)
  