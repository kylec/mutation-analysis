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
  universe = c()
  venn_img = paste("venn", i,"tiff", sep=".")
  cols<-c("cornflowerblue", "green", "yellow", "darkorchid1", "red")
  heatmap_img = paste("heatmap", i, "png", sep=".")
  
  # add genes of samples to a list
  for (i in 1:length(samples)) {
    file_name = strsplit(sample[i], "_")[[1]][1]
    fin = paste(file_name,"ann.tsv",sep=".")
    a = read.table(fin, sep="\t", header=T, stringsAsFactors=F)
    #b = a[grep("exonic", a$region_type),]
    b = a[grep("nonsyn", a$mut_type),]
    # max sp
    #max_sp = max(na.omit(b$SP))
    genes[[i]] = b$genename
    universe = append(universe, b$genename)
  }
  universe = sort(unique(universe))
  
  # create mutation heatmap matrix
  # Generate a matrix, with the sets in columns and possible letters on rows
  Counts <- matrix(0, nrow=length(universe), ncol=length(samples))
  
  # Populate the said matrix
  for (i in 1:length(universe)) {
    for (j in 1:length(samples)) {
      Counts[i,j] <- universe[i] %in% genes[[j]]
    }
  }
  
  # Name the columns with the sample names
  colnames(Counts) <- samples
  rownames(Counts) = universe
  
  # filter to non-singleton hits
  format_data = melt(Counts[rowSums(Counts==1)>1,])
  
  # plot mutation heatmap
  names(format_data) = c("variable", "sample", "value")
  base_size <- 12
  
  (p <- ggplot(format_data, aes(sample, variable)) + geom_tile(aes(fill = value),colour = "black") 
   + scale_fill_gradient(low = "white",high = "steelblue") 
   + theme_grey(base_size = base_size) 
   + labs(x = "", y= "") + scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) 
   + theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = base_size , angle = 90, hjust = 0, colour = "grey50")))  
  
  ggsave(heatmap_img, p)
  
  # plot venn
  library(VennDiagram)
  venn.diagram(genes,  category.names = samples,
               filename=venn_img, col = "transparent", fill = cols[1:length(samples)], alpha = 0.50, scaled=TRUE,
               cat.cex=.9, cat.dist = rep(0.1, length(samples)))
  
  if (1 != 1) {
  #### TEST ####
  ####  gplot
  Reduce(intersect, genes)
  library(gplots)
  c=venn(list(Vilar11=genes[[1]],Vilar12=genes[[2]], Vilar22=genes[[3]], Vilar23=genes[[4]], Vilar24=genes[[5]]))
  
  #### limma
  vennDiagram(vennCounts(Counts), circle.col=cols)
  }
  
} # end of sample groupings loop
