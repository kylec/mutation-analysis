library(ggplot2)
library(reshape2)

# get a list of genes for each samples
# create a whole set of genes

# read list of sample groupings
f = system("cat /Users/kyle_air/Dropbox/lab_vilar/fap/samples_group.txt", intern=T)

for (i in 1:length(f)) {
  
  # get a list of samples (Vilar02_N) of a patient
  samples = strsplit(f[i], "\t")[[1]]
  genes = list() # a list of samples , which has a list of genes
  venn_img = paste("venn", i,"tiff", sep=".")
  cols<-c("cornflowerblue", "green", "yellow", "darkorchid1", "red")
 
  
  # add genes of samples to a list
  for (i in 1:length(samples)) {
    file_name = strsplit(samples[i], "_")[[1]][1]
    fin = paste(file_name,"ann.tsv",sep=".")
    a = read.table(fin, sep="\t", header=T, stringsAsFactors=F)
    #b = a[grep("exonic", a$region_type),]
    b = a[grep("nonsyn", a$mut_type),]
    # max sp
    #max_sp = max(na.omit(b$SP))
    genes[[i]] = b$genename
  }
  
  
  # plot venn
  library(VennDiagram)
  venn.diagram(genes,  category.names = samples,
               filename=venn_img, col = "transparent", fill = cols[1:length(samples)], alpha = 0.50, scaled=TRUE,
               cat.cex=.9, cat.dist = rep(0.15, length(samples)))
  
  #### TEST ####
  ####  gplot
  #3Reduce(intersect, genes)
  #library(gplots)
  #c=venn(list(Vilar11=genes[[1]],Vilar12=genes[[2]], Vilar22=genes[[3]], Vilar23=genes[[4]], Vilar24=genes[[5]]))
  #### limma
  #vennDiagram(vennCounts(Counts), circle.col=cols)

  
} # end of sample groupings loop
