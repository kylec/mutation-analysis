library(cummeRbund)
#library(graphics)
library(amap)
library(gplots)
clust_method = "pearson"

plot_heatmap = function(path, heatmap_name) {
  setwd(path)
  cat("reading cuffdata...")
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
  # sigGeneIds <- getSig(t_vs_n_cuff_data, level="genes", alpha=0.005)
  
  #cuffdiff object
  cat("creating cuffdiff object...")
  t_vs_n_diff_genes <- diffData(genes(t_vs_n_cuff_data))
  # TODO: significant gene id - why not use getSig? getSig return less genes
  t_vs_n_diff_genes_sig_ids <- t_vs_n_diff_genes[t_vs_n_diff_genes$q_value <= 0.05 & abs(t_vs_n_diff_genes$log2_fold_change) >= 1,]$gene_id
  t_vs_n_de_genes <- getGenes(t_vs_n_cuff_data, t_vs_n_diff_genes_sig_ids)
  #csHeatmap(t_vs_n_de_genes,cluster='both', replicates=T)
  
  t_matrix <- repCountMatrix(t_vs_n_de_genes)
  
  #TODO: check if up/down regulated genes make sense because they foldchange is normal/tumor
  clust.genes<-hcluster(x=as.matrix(t_matrix), method="correlation", link="average")
  clust.samples<-hcluster(x=t(as.matrix(t_matrix)), method="euclidean", link="average")
  STATUS =c(rep("white",length(grep("normal",clust.samples$labels))), rep("black", length(grep("polyp",clust.samples$labels))))
  png(paste0(path,heatmap_name), width=1000, height=800, res=100)
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
path="~/Analysis/fap/rna_seq/cdout/colon_polyp,colon_normal/"
heatmap_name="colon_polyp_normal_heatmap.png"
plot_heatmap(path, heatmap_name)

### duodenum polyp
path="~/Analysis/fap/rna_seq/cdout/duodenum_polyp,duodenum_normal/"
heatmap_name="duodenum_polyp_normal_heatmap.png"
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


