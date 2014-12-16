# read list of sample groupings
#library(ggplot2)
f = system("cat ~/Dropbox/lab_vilar/fap/samples_group.txt", intern=T)

for (i in 1:length(f)) {
  
  # get a list of samples (Vilar02_N) of a patient
  samples = strsplit(f[i], "\t")[[1]]
  genes = list() # a list of samples , which has a list of genes
  cols<-c("cornflowerblue", "green", "yellow", "darkorchid1", "red")
  clone_img = paste("clone", i, "png", sep=".")
  af_img = paste("af", i, "png", sep=".")
  # af plot list holder for ggplot
  plot_list = c() 
  
  # plot a list of samples per patient
  png(clone_img, width=800, height=400)
  par(mfrow=c(1,length(samples)))
  for (i in 1:length(samples)) {
    file_name = strsplit(samples[i], "_")[[1]][1]
    fin = paste(file_name,"sps.ann.tsv",sep=".")
    a = read.table(fin, sep="\t", header=T, stringsAsFactors=F)
    
    sp = sort(unique(a$SP[!is.na(a$SP)]))
    
    plot(seq(1:length(sp)), sp, main=samples[i], xlab="clones #", ylab="Tumor %", cex=2, pch=16, ylim=c(0,1) )
   
    # AF plots
    b = a[which(a$PN_B==0 & a$f!="NA" & a$f <= 1),]
    #b=b[grep("exonic", b$region_type),]
    p = ggplot(b) +
      geom_density(aes(x=AF_Tumor, fill="AF"), alpha=.3) + 
      geom_density(aes(x=f, fill="Adj. AF"), alpha=.3) + 
      #geom_histogram(aes(x=AF_Tumor), binwidth=.01, fill="blue", alpha=.2) +
      #geom_histogram(aes(x=f), binwidth=.01, fill="red", alpha=.2) +
      xlab("Allele_Fraction") + ggtitle(samples[i]) + theme(text=element_text(size=18))
    plot_list[[i]] = p
  }
  dev.off()
  library(grid)
  library(gridExtra)
  names(plot_list) = samples
  png(af_img, width=2000, height=1000)
  arg_list <- c(plot_list, list(nrow=1, ncol=length(samples)))
  do.call(grid.arrange, arg_list)
  dev.off()
}