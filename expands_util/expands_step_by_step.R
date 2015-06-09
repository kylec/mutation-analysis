library(expands)

args <- commandArgs(TRUE)

basename = args[1]
snvFile = paste(basename, "combined", sep=".")
cnFile = paste(basename, "copynumber.dnacopy", sep=".")

# read snv file
cat("reading snv.")
snv = read.table(snvFile, sep="\t", header=T)
snv = data.matrix(snv)

# read cnv file
cat("reading cnv.")
cbs = read.table(cnFile, sep="\t", header=F)
cbs = cbs[,c(1,2,3,5)]
colnames(cbs) = c("chr", "startpos","endpos", "CN_Estimate")
head(cbs)
# remove log2 ratio
cbs$CN_Estimate = (2^cbs$CN_Estimate)*2
cbs = data.matrix(cbs)

# step by step expands
cat("running expands.")
# assign snp to cnv segments
dm=assignQuantityToMutation(snv,cbs,"CN_Estimate")
dm = dm[which(dm[,"CN_Estimate"]!="NA"),]
max_PM=6; maxScore=2.5; precision=0.018;
plotF=1
##the name of the sample
snvF=basename
cfd=computeCellFrequencyDistributions(dm, max_PM, precision)
toUseIdx=which(apply(is.finite(cfd$densities),1,all) )
SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], precision, label=snvF)
aM= assignMutations( dm, SPs,cfd$densities)
png(paste(basename, "png", sep="."))
plotSPs(aM$dm, snvF,cex=1)
dev.off()
aQ=assignQuantityToSP(cbs, aM$dm)
tr=buildPhylo(aQ,snvF)
plot(tr,cex=2.5)
write.table(aM$dm, file=paste(snvF,"sps", sep="."))