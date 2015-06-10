library(knitr)
# plot cnv

#files=list.files(path=".", pattern="called.dnacopy$")
#files=list.files(path=".", pattern="called.dnacopy.filtered$")
#files=list.files(path=".", pattern="copynumber.dnacopy$")
#files=list.files(path=".", pattern="called.dnacopy.merged$")

# set chr length, range , color
chrLen = data.frame(chr=seq(1,24,by=1), length=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,155270560,59373566) )
sum=0
chrNum=length(chrLen$length)
padding=c(0)
for (i in 2:chrNum) {
  sum = sum + chrLen$length[i-1]
  padding = append(padding, sum)
}
colors=rep("black", chrNum)
colors[seq(1,chrNum, by=2)] = "red" 
chrLen = cbind(chrLen, padding, colors)

# plot cnv
read_data = function(files) {
  mean_lrr = NULL

  for (file in files) {
    # merge paddings and color to seg file
    name = file
    a = read.table(file, sep="\t", header=T)
    colnames(a) = c('chr', 'start','end','seglen','lrr')
    # remove chr
    a$chr = as.numeric(gsub("chr","", a$chr))
    b = merge(a, chrLen, by=c('chr'))
    
    padstart = b$start+b$padding
    padend = b$end+b$padding
    b = cbind(b, padstart, padend)
    
    # plot cnv 
    png(paste0(name,'.png'), height=400, width=1200)
    plot(c(b$padstart, b$padend), c(b$lrr, b$lrr), col=b$colors, cex=0.2, main=name, ylab="log2ratio", ylim=c(-2,2))
    par(new=T)
    abline(v=b$padding)
    dev.off()
    
  
    # mean lrr 
    mean_lrr = append(mean_lrr, mean(b$lrr))
  }
  return(mean_lrr)
}

# main 
setwd("/Users/kchang3/Analysis/fap/varscan")
files=list.files(path=".", pattern="called.dnacopy$")
mean_lrr = read_data(files)

setwd("/Users/kchang3/Analysis/fap/varscan-adj")
files=list.files(path=".", pattern="called.dnacopy$")
mean_lrr_adj = read_data(files)

setwd("/Users/kchang3/Analysis/fap/varscan-adj-tgt")
files=list.files(path=".", pattern="called.dnacopy$")
mean_lrr_adj = read_data(files)

cnv_summary = data.frame(sample=substr(files, 1, 7), mean_lrr=mean_lrr, mean_lrr_adj=mean_lrr_adj)

setwd("/Users/kchang3/Analysis/tcga_coadread/varscan")
files=list.files(path=".", pattern="called.dnacopy$")
mean_lrr_adj = read_data(files)

setwd("/Users/kchang3/Analysis/tcga_coadread/varscan-adj-tgt")
files=list.files(path=".", pattern="called.dnacopy$")
mean_lrr_adj = read_data(files)

setwd("/Users/kchang3/Analysis/tcga_coadread/cnv")
files=list.files(path=".", pattern=".seg$")
mean_lrr_adj = read_data(files)



