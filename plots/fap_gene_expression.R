library(cummeRbund)
#library(graphics)
library(amap)
library(gplots)
clust_method = "pearson"
# TODO: significant gene id - why not use getSig? getSig return less genes
# TODO: check if up/down regulated genes make sense because they foldchange is normal/tumor
# TODO: we are counting each sample as a replicate of either polyp  or normal , is it right?

plot_heatmap = function(path, heatmap_name) {
  setwd(path)
  t_vs_n_cuff_data <- readCufflinks()
  # runInfo(t_vs_n_cuff_data)
  # dispersionPlot(genes(t_vs_n_cuff_data))
  # csDensity(genes(t_vs_n_cuff_data))
  # csScatter(genes(t_vs_n_cuff_data), 'colon_polyp', 'colon_normal')
  # csDendro(genes(t_vs_n_cuff_data), replicates=TRUE)
  # replicates(t_vs_n_cuff_data) # describes replicates (can get sample names with this command)
  # csBoxplot(genes(t_vs_n_cuff_data), replicates=TRUE)
  # csScatterMatrix(genes(t_vs_n_cuff_data))
  # csVolcano(genes(t_vs_n_cuff_data),'colon_polyp','colon_normal', alpha=0.05, showSignificant=T)
  
  # http://seqanswers.com/forums/showthread.php?t=19278
  # stringent filter for sig de genes id
  sigGeneIds <- getSig(t_vs_n_cuff_data, level="genes", alpha=0.05)
  
  # return gene object of a cuffset, then return diff gene dataframe - value1 and value2 are average FPKM in genes.read_group_tracking
  t_vs_n_diff_genes <- diffData(genes(t_vs_n_cuff_data),features=TRUE)
  # anthony filter for sig de genes id
  t_vs_n_diff_genes_sig_ids <- t_vs_n_diff_genes[t_vs_n_diff_genes$q_value <= 0.05 & abs(t_vs_n_diff_genes$log2_fold_change) >= 1,]$gene_id
  
  # write de gene file by 1) getSig, 2) anthony's filter
  write.table(t_vs_n_diff_genes[ t_vs_n_diff_genes$gene_id %in% sigGeneIds, ], file="de_gene_getSig.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  write.table(t_vs_n_diff_genes[ t_vs_n_diff_genes$gene_id %in% t_vs_n_diff_genes_sig_ids, ], file="de_gene_anthony.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  
  # returns "internal_scaled_frags" in file genes.read_group_tracking - Estimated number of fragments originating from the object, after transforming to the internal common count scale (for comparison between replicates of this condition.)
  # cuffgeneset using anthony's filter for de genes
  t_vs_n_de_genes <- getGenes(t_vs_n_cuff_data, t_vs_n_diff_genes_sig_ids)
  #csHeatmap(t_vs_n_de_genes,cluster='both', replicates=T)
  t_matrix <- repCountMatrix(t_vs_n_de_genes)
 
  # cluster fragments count
  clust.genes<-hcluster(x=as.matrix(t_matrix), method="correlation", link="average")
  clust.samples<-hcluster(x=t(as.matrix(t_matrix)), method="euclidean", link="average")
  STATUS =c(rep("white",length(grep("normal",clust.samples$labels))), rep("black", length(grep("polyp",clust.samples$labels))))
  png(paste0(path,heatmap_name,"_anthony.png"), width=1000, height=800, res=100)
  heatmap.2(x=as.matrix(t_matrix),
            Rowv=as.dendrogram(clust.genes),
            Colv=as.dendrogram(clust.samples),
            mar=c(15,1),
            col=redblue(75),
            key=TRUE,
            # symkey=FALSE,
            # density.info="none",
            trace="none",
            sub="",
            ColSideColors=STATUS,
            scale="row",
            #main=paste(output_prefix, "differentially expressed genes", sep=" ")
            # sepwidth=c(0.00001,0.00001),
            # sepcolor="#DDDDDD",
            # colsep=1:ncol(esetSel),
            # rowsep=1:nrow(esetSel),
  )
  dev.off()
  
  # getSig genes clustering using stringent alpha cutoff
  t_vs_n_de_genes <- getGenes(t_vs_n_cuff_data, sigGeneIds)
  t_matrix <- repCountMatrix(t_vs_n_de_genes)
  # cluster fragments count
  clust.genes<-hcluster(x=as.matrix(t_matrix), method="correlation", link="average")
  clust.samples<-hcluster(x=t(as.matrix(t_matrix)), method="euclidean", link="average")
  STATUS =c(rep("white",length(grep("normal",clust.samples$labels))), rep("black", length(grep("polyp",clust.samples$labels))))
  png(paste0(path,heatmap_name,"_getSig.png"), width=1000, height=800, res=100)
  heatmap.2(x=as.matrix(t_matrix),
            Rowv=as.dendrogram(clust.genes),
            Colv=as.dendrogram(clust.samples),
            mar=c(15,1),
            col=redblue(75),
            key=TRUE,
            # symkey=FALSE,
            # density.info="none",
            trace="none",
            sub="",
            ColSideColors=STATUS,
            scale="row",
            #main=paste(output_prefix, "differentially expressed genes", sep=" ")
            # sepwidth=c(0.00001,0.00001),
            # sepcolor="#DDDDDD",
            # colsep=1:ncol(esetSel),
            # rowsep=1:nrow(esetSel),
  )
  dev.off()
}

### colon_polyp
path="~/Analysis/fap/rnaseq/cdout/colon_polyp,colon_normal/"
heatmap_name="colon_polyp_normal_heatmap"
plot_heatmap(path, heatmap_name)

### duodenum polyp
path="~/Analysis/fap/rnaseq/cdout/duodenum_polyp,duodenum_normal/"
heatmap_name="duodenum_polyp_normal_heatmap"
plot_heatmap(path, heatmap_name)

### lynch
path="~/Analysis/lynch/rnaseq-human/cdout/colon_polyp,colon_normal/"
heatmap_name="colon_polyp_normal_heatmap"
plot_heatmap(path, heatmap_name)

### polyp vs normal
#TODO: filter polyp matrix for colon and duodenum gene signature

###  gene converter
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


