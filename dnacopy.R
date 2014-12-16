library(DNAcopy)

# dnacopy.R input output

args <- commandArgs(TRUE)
cn <- read.table(args[1],header=T)
CNA.object <-CNA( genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file=args[2], row.names=F, col.names=F, quote=F, sep="\t")
# write p values for segment to be used in varscan mergeSegments.pl
pseg = segments.p(segs, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)
write.table(pseg, file=args[3], row.names=F, col.names=T, quote=F, sep="\t")
