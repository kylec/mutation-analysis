# convert mouse gene to human
# usage: Rscript convertGene.R [gene_list.txt(no header)] [human/mouse]
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
library(plyr)

args <- commandArgs(TRUE)

diff_gene = args[1]
gene_type = args[2]

#test
#setwd("/Users/kchang3/Analysis/tsai_skin/skin_mice")
#diff_gene = "/Users/kchang3/Analysis/tsai_skin/skin_mice/diff_gene.txt"
#gene_type = "mouse"
#genes = c("AW112010","AW209491","AW551984","Aaas","A930018M24Rik","Aar2","Sf3b6")
#genes = c("DSTN","TP53","PTEN","BRCA1")

# read gene list
df = read.table(diff_gene, header=T, sep="\t")
df$gene = as.character(df$gene)
#TODO: modify for other gene column
genes = df$gene

# load ensembl genes
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# human to mouse genes
if (gene_type == "human") {
  genes_df = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genes ,mart = human, attributesL = c("mgi_symbol"), martL = mouse)
  colnames(genes_df) = c("gene","newgene")
  write.table(genes_df, "gene_lookup.txt", col.names=F, row.names=F, sep="\t", quote=F)
  # mouse to human genes
} else {
  genes_df = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes ,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  colnames(genes_df) = c("gene","newgene")
  write.table(genes_df, "gene_lookup.txt", col.names=F, row.names=F, sep="\t", quote=F)
}

# input genes without human id (not in biomart lookup)
write.table(genes[! genes %in% genes_df$gene], "not_found_genes.txt", col.names=F, row.names=F, sep="\t", quote=F)

# find mouse gene map to multiple human gene
dup_genes = as.character(count(genes_df$gene)[count(genes_df$gene)$freq>1,]$x)

# get a data frame of biomrt genes lookup output
dup_genes_df = genes_df[genes_df$gene %in% dup_genes,]

#TODO: hack to mouse genes if uppercase mouse name == human name , otherwise I can't resolve a mouse gene that has multiple unrelated human genes names
# which means i'm ignoring some mouse genes
dup_genes_df$ucgene = toupper(dup_genes_df$gene)

# dedup genes
dedup_genes_df = dup_genes_df[dup_genes_df$newgene == dup_genes_df$ucgene, c("gene","newgene")]
# save dup genes that can't be resolved
write.table(dup_genes_df[! dup_genes_df$gene %in% dedup_genes_df$gene, c("gene","newgene")], "dup_genes.txt", col.names=F, row.names=F, sep="\t", quote=F)

# combine uniquely mapped genes and de-dup genes 
dedup_genes_df = rbind(genes_df[! genes_df$gene %in% dup_genes,], dedup_genes_df)
#write.table(dedup_genes_df, "test.txt", col.names=F, row.names=F, sep="\t", quote=F)

# write diff gene with mouse and human gene name
#dim(df[df$gene %in% dedup_genes_df$gene, ])
dff = merge(df, dedup_genes_df, by=c("gene"))
dff$mouse_gene = dff$gene
dff$gene = dff$newgene 
dff$newgene = NULL
write.table(dff, "human_diff_gene.txt", row.names=F, sep="\t", quote=F)

# head(df)
# head(genes)
# dim(genes)
# dim(dff)
# dim(dedup_genes_df)
# length(unique(dedup_genes_df$gene))
# length(unique(dup_genes_df$gene))
# length(unique(genes_df$gene))
# 
# intersect(unique(dedup_genes_df$gene), unique(dup_genes_df$gene))

#TODO: these genes appeared AFTER biomart lookup, which were not in original gene input . I have no idea why
setdiff(dedup_genes_df$gene, dff$gene) # "Rpl24"    "Tmem185b" "Kctd12" 
setdiff(dff$gene, dedup_genes_df$gene)

# library(biomaRt)
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# genes = c("Zfp286", "Tmx2")
# genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes ,mart = mouse, attributesL = c("hgnc_symbol","chromosome_name", "start_position"), martL = human, uniqueRows=T)
