# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library("DESeq2")

setwd("/Users/kchang3/Analysis/fap/rnaseq-human/htseq")

# can merge individual sample files (i.e. ctrl1.counts, ctrl2.counts, etc.)
sampleFiles = system("ls tophat*.txt", intern=T)
samples = read.table("sample_groups/samples.txt", header=T, sep="\t", stringsAsFactors = F)
samples$sample_name = gsub("_","",samples$sample_name)
samples$patientid = paste0(samples$patient_id,"_",samples$sample_name,"_",samples$Type,"_",samples$Localization)

head(samples)

# view sampleFiles
sampleFiles
index = gsub("_","",substr(sampleFiles,15,19))

# can designate different batches of samples (i.e. different sequencers,
# PE vs SE, different library preps (eg. BATCH1 vs BATCH2))
#sampleBatch <- c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1", "Batch2","Batch2","Batch2","Batch2")

samples = samples[match(index, samples$sample_name),]

condition = samples$Type
# sampleTable <- data.frame(sampleName = sampleFiles,
#                           fileName = sampleFiles,
#                           condition = sampleCondition,
#                           Batch = sampleBatch)
sampleTable <- data.frame(sampleName = samples$patientid, fileName = sampleFiles, condition = condition)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable, directory = ".", design = ~condition) 
colData(ddsHTseq)$condition<-factor(colData(ddsHTseq)$condition, levels=unique(condition))

dds <- DESeq(ddsHTseq)

res <- results(dds)
# order results by padj value (most significant to least)
res <- res[order(res$padj),]
head(res)
res[row.names(res)=="CNOT3",]

plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld)),file='DATE-DESeq2-rlog-transformed-counts.txt', sep='\t')
write.table(as.data.frame(assay(vsd)),file='DATE-DESeq2-vst-transformed-counts.txt', sep='\t')
                        
# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,"DATE-DESeq2_variance_stabilizing.png")
dev.off()

library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png("DESeq2_heatmap1.png", height=800, width=800)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6))
dev.off()
png("DESeq2_heatmap2.png", height=800, width=800)
heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
dev.off()
png("DESeq2_heatmap3.png", height=800, width=800)
heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition,sampleFiles , sep=" : "))

hc <- hclust(distsRL)
png("deseq2_heatmaps_samplebysample.png", height=800, width=800)
heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()


t_matrix = assay(rld)

# cluster fragments count
# clust.genes<-hcluster(x=as.matrix(t_matrix), method="correlation", link="average")
# clust.samples<-hcluster(x=t(as.matrix(t_matrix)), method="euclidean", link="average")
# STATUS =c(rep("white",length(grep("normal",clust.samples$labels))), rep("black", length(grep("polyp",clust.samples$labels))))
# png(paste0(path,heatmap_name,"_anthony.png"), width=1000, height=800, res=100)
# heatmap.2(x=as.matrix(t_matrix),
#           Rowv=as.dendrogram(clust.genes),
#           Colv=as.dendrogram(clust.samples),
#           mar=c(15,1),
#           col=redblue(75),
#           key=TRUE,
#           # symkey=FALSE,
#           # density.info="none",
#           trace="none",
#           sub="",
#           ColSideColors=STATUS,
#           scale="row",
#           #main=paste(output_prefix, "differentially expressed genes", sep=" ")
#           # sepwidth=c(0.00001,0.00001),
#           # sepcolor="#DDDDDD",
#           # colsep=1:ncol(esetSel),
#           # rowsep=1:nrow(esetSel),
# )
