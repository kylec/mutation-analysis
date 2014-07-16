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

# Filter loh > .5 in snv file
loh = snv[which(snv$PN_B==1 & snv$AF_Tumor > .5),]
# max sampling number of loh markers
sampleNums = floor(seq(1, dim(loh)[1], length.out=5))
outputName=basename

# loop through sampling number
counter = 0
for (sampleNum in sampleNums) {
counter = counter + 1
# append basename with the sampling number
outputName = paste(basename, counter, sep="-" )

if (dim(loh)[1] > 1) {  
  print(paste("sampling =", sampleNum, sep=" "))
  
  # original
  # keep filtered loh (x sd from mean) and somatic
  #mean_af = mean(loh$AF_Tumor)
  #sd_af = sd(loh$AF_Tumor)
  #snv = snv[which(snv$PN_B==1 & (snv$AF_Tumor <= mean_af-sd_num*sd_af | snv$AF_Tumor >= mean_af+sd_num*sd_af) | snv$PN_B ==0), ]
  
  # sample a subset of the loh markers
  sample_loh = loh[sample(nrow(loh), sampleNum),]
  sample_loh$AF_Tumor = mean(sample_loh$AF_Tumor)
  
  # combined snv + processed loh markers
  combined = rbind(snv[snv$PN_B==0,], sample_loh)
  
}
print(paste("size =", dim(combined), sep=" "))

if (FALSE) {
combined = data.matrix(combined)
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
runExPANdS(combined, cbs, maxScore=2.5, max_PM=6, precision=NA,plotF=2,snvF=outputName,maxN=8000,region=NA)
}
} # end of sampling