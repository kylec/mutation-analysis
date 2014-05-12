library(expands)

args <- commandArgs(TRUE)

# inputs
basename = args[1]
#basename = "Vilar05"
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
# remove log2 ratio
cbs$CN_Estimate = 2*(2^cbs$CN_Estimate)
cbs = data.matrix(cbs)

# run expands
cat("running expands.")
runExPANdS(snv, cbs, maxScore=2.5, max_PM=6, precision=NA,plotF=2,snvF=basename,maxN=8000,region=NA)
