library(CHAT)

## functions
ParseVCF <- function(filename,tumor,normal,AD,thr.cov){
  vcf<-read.table(filename,sep='\t',header=F,stringsAsFactors=F)
  # remove multi-allelic positions 
  vcf = vcf[which(nchar(vcf$V4) == 1 & nchar(vcf$V5) == 1), ]
  sampleid=sub('.raw.vcf','',filename)
  
  # determine tumor and normal column 
  header = unlist(strsplit(system(paste("grep CHROM", filename), intern=T), "\t"))
  
  # special - not every sample is Vilar## , just grep by number
  #tumor = grep(tumor, header)
  #normal = grep(normal, header)
  tumor = grep(gsub("Vilar", "", tumor), header)
  normal = grep(gsub("Vilar", "", normal), header)

  # index for columns that have allele data
  vv.a=which(vcf[,tumor]!='./.'&vcf[,normal]!='./.')
  vcf=vcf[vv.a,]
  
  # returns AP field  GT:AD:DP:GQ:PL  1/1:0,8:8:24:291,24,0
  # each row has 2 element, ref count and alt count
  y=unlist(lapply(strsplit(vcf[,tumor],':'),function(x)x[[AD]]))
  y=strsplit(y,',')
  y0=unlist(lapply(strsplit(vcf[,normal],':'),function(x)x[[AD]]))
  y0=strsplit(y0,',')
  
  # pick rows that have ref and alt count?
  vv.bi=which(unlist(lapply(y,length))==2&unlist(lapply(y0,length))==2)
  y=y[vv.bi]
  y0=y0[vv.bi]
  vcf=vcf[vv.bi,]
  
  # a = ref count, b = alt count, s = total depth
  a=as.numeric(unlist(lapply(y,function(x)x[[1]])))
  b=as.numeric(unlist(lapply(y,function(x)x[[2]])))
  s=a+b
  a0=as.numeric(unlist(lapply(y0,function(x)x[[1]])))
  b0=as.numeric(unlist(lapply(y0,function(x)x[[2]])))
  s0=a0+b0
  #gt.n=unlist(lapply(strsplit(vcf[,normal],':'),function(x)x[[1]]))
  #gt.t=unlist(lapply(strsplit(vcf[,tumor],':'),function(x)x[[1]]))
  # make sure the normal ref and alt counts are > 0
  # normal and tumor total coverage > threshold
  vv.germ=which(a0>0&b0>0)
  vv.cov=which(s>=thr.cov&s0>=thr.cov)
  vv=intersect(vv.germ,vv.cov)
  
  seg.mat=cbind(sampleid,gsub("chr","",vcf[vv,1]),vcf[vv,2],log2(s/s0)[vv],b[vv]/s[vv],b0[vv]/s0[vv])
  colnames(seg.mat)=c('sampleID','chr','pos','LRR','BAF','BAF-n')
  return(seg.mat)
}

getSegChr.Seq <- function(seg.mat,bin=1000,cbs=TRUE,thr.hets=0.15){
  sampleid=sub('.raw.vcf','',filename)
  dd.dat=c()
  id=seg.mat[1,1]
  for(cc in 1:22){
    vv=which(seg.mat[,2]==cc)
    bb.chr=seg.mat[vv,c(2,2,3,5,6)]
    ll.chr=cbind(seg.mat[vv,c(2,2,3,4)],rep(0,length(vv)))
    mode(bb.chr)=mode(ll.chr)='numeric'
    colnames(bb.chr)[4:5]=colnames(ll.chr)[4:5]=c(sampleid,paste(sampleid,'-normal',sep=''))
    if (cbs) { 
      dat.chr=getSegChr.CBS(bb.chr,ll.chr,sam.col=4,thr.hets=thr.hets,data.type='log')
    } else { 
      dat.chr=getSegChr(bb.chr,ll.chr,sam.col=4,thr.hets=thr.hets,data.type='log')
    }
    dd.dat=rbind(dd.dat,dat.chr)
  }
  return(dd.dat)
}

## main
project_dir="~/Projects/fap/"
setwd(paste0(project_dir, "chat"))

dd.dat = NULL
samples = read.table(paste0(project_dir, "pairs.txt"),header=F,stringsAsFactors=F)
# loop through patients.txt
for (i in 1:length(samples$V2)) {
  tumor=samples[i,]$V2
  normal=samples[i,]$V3
  AD=2
  thr.cov=20
  filename = paste0(tumor, ".raw.vcf")
  seg.mat = ParseVCF(filename, tumor, normal, AD, thr.cov)

  dd.tmp = getSegChr.Seq(seg.mat)
  dd.dat = rbind(dd.dat, dd.tmp)

}


#save segmentation data as .rdata
sampleid="fap"
seg.data=paste0(sampleid,".Rdata")
seg_sAGP.data=paste0(sampleid,"_sAGP.Rdata")
agp.txt = paste0(sampleid, "_AGP.txt")

rownames(dd.dat) = dd.dat[,"id"]
dd.dat = dd.dat[,c(2:8)]
mode(dd.dat)='numeric'
save(dd.dat, file=seg.data)

#AGP inference
para <- getPara()
para$datafile <- seg.data 
para$thr.penalty <- 300
para$savefile <- agp.txt
getAGP(para=para)

#sAGP
para.s <- getPara.sAGP()
para.s$inputdata <- seg.data
para.s$purityfile <- agp.txt
para.s$savedata <- seg_sAGP.data
getsAGP(para=para.s)