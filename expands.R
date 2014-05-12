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
# Filter loh in snv file
head(snv)
b = snv[which(snv$PN_B==1),]
sd_num = 1.5
if (dim(b)[1] > 1) {
  mean_af = mean(b$AF_Tumor)
  sd_af = sd(b$AF_Tumor)
  print(paste(basename, round(mean_af,3), round(sd_af,3), sep="\t"))
  # keep filtered loh (x sd from mean) and somatic
  snv = snv[which(snv$PN_B==1 & (snv$AF_Tumor <= mean_af-sd_num*sd_af | snv$AF_Tumor >= mean_af+sd_num*sd_af) | snv$PN_B ==0), ]
}

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
