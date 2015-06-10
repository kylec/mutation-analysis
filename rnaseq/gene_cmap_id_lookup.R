
# compare gene list
dir="/Users/kchang3/Analysis/fap/rnaseq-human"
setwd(dir)

file1 = paste0(dir,"/cdout/colon_normal_polyp/up.txt")
up1 = read.table(file1)

file2 = paste0(dir,"/anthony/colon_normal_polyp/up.txt")
up2 = read.table(file2)

dim(up1)
dim(up2)

length(intersect(up1$V1, up2$V1))
length(setdiff(up1$V1, up2$V1))
length(setdiff(up2$V1, up1$V1))


# convert gene symbol to affymetrix id 
#source("http://bioconductor.org/biocLite.R")
#biocLite("hgu133a.db")
library(hgu133a.db)
library(annotate)
x <- hgu133aSYMBOL
head(x)
# Get the probe identifiers - gene symbol mappings
mapped_probes <- mappedkeys(x)
# Convert to a dataframe
genesym.probeid <- as.data.frame(x[mapped_probes])
head(genesym.probeid)
write.table(genesym.probeid, "array_gene_id.txt", col.names=F, row.names=F, quote=F, sep="\t")
#read up/down list
setwd("/Users/kchang3/Analysis/fap/rnaseq-human/anthony/colon_normal_polyp")
up= read.table("up.txt")
colnames(up) = "symbol"
d = merge(up, genesym.probeid, by=c("symbol") )
write.table(unique(d$probe_id), "up_array_id.grp", col.names=F, row.names=F, quote=F, sep="\t")
down= read.table("down.txt")
colnames(down) = "symbol"
d = merge(down, genesym.probeid, by=c("symbol") )
write.table(unique(d$probe_id), "down_array_id.grp", col.names=F, row.names=F, quote=F, sep="\t")

