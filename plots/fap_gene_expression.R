library(cummeRbund)

setwd("~/projects/fap/rna_seq")
cuff<-readCufflinks()

# add features
annot<-read.table("gene_annotation.tab",sep="\t",header=T,na.string="-") 
addFeatures(cuff,annot,level="genes")

disp<-dispersionPlot(genes(cuff))
genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
dens<-csDensity(genes(cuff))
dens

# sig gene set level
mySigGeneIds = getSig(cuff,alpha=0.05,level='genes')
myGenes<-getGenes(cuff,mySigGeneIds)
h<-csHeatmap(myGenes,cluster='both')
h

hh = getSigTable(cuff,alpha=0.05,level='genes')
hh

gene.diff<-diffData(genes(cuff))
dim(gene.diff[gene.diff$significant=='yes',])

## search db
# load gene exp 
# load gene features (fpkm(myGenes,features=T)
# merge gene exp and gene name in gene features
# if mouse, merge human names to mouse names

# load mouse gene exp
# replace mouse gene name with human gene



source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
human_genes = c("DSTN","TP53","PTEN","BRCA1")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
getLDS(attributes = c("hgnc_symbol"),
       filters = "hgnc_symbol", values = human_genes ,mart = human,
       attributesL = c("mgi_symbol"), martL = mouse)


