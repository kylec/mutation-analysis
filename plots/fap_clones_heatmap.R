# list of all mutated genes

library(ggplot2)
library(reshape2)

# read list of sample groupings
f = system("cat /Users/kyle_air/Dropbox/lab_vilar/fap/samples_group.txt", intern=T)

# all the genes mutated across all samples
crc_genes = system("cat crc_genes.txt", intern=T)
#universe = mutated_genes
#heatmap_img = paste("fap_clones_mut_heatmap", "png", sep=".")

#master_heatmap = NULL
#xaxis_text_colors = NULL
#color_counter = 0

# get a list of genes for each samples
# create a whole set of genes
# each patient, each sample
for (i in 1:length(f)) {
  
  #### get universe genes by patient
  master_heatmap = NULL
  xaxis_text_colors = NULL
  yaxis_text_colors= NULL
  color_counter = 0
  
  s = gsub("\t", "*.tsv ", f[i])
  s = gsub("_[A-Z]", "", s)
  cmd = paste("awk -F \"\t\" \'$4==0\' ", s, "*.tsv | cut -f13-15 | grep exonic | egrep \"stop|nonsyn\" | cut -f2 | sort -u", sep="")
  mutated_genes = system(cmd, intern=T)
  universe = mutated_genes
  heatmap_img = paste(i, "heatmap", "png", sep=".")
  
  # get a list of samples (Vilar02_N) of a patient
  samples = strsplit(f[i], "\t")[[1]]
  
  # add genes of samples to a list
  for (j in 1:length(samples)) {
    file_name = strsplit(samples[j], "_")[[1]][1]
    fin = paste(file_name,"ann.tsv",sep=".")
    sample_data = read.table(fin, sep="\t", header=T, stringsAsFactors=F)
    
    # each sample's sp, get genes
    sps = sort(unique(na.omit(sample_data$SP)))
    
    names=c()
    for (q in 1:length(sps)) {
      names = append(names, paste(samples[j], round(sps[q],2), sep="_"))
    }
    
    # heat map matrix with tumor % 
    sample_heatmap <- matrix(0, nrow=length(universe), ncol=length(sps))
    colnames(sample_heatmap) <- names
    rownames(sample_heatmap) = universe
    
    # alternate x axis label color to separate polyp groups
    if (color_counter %% 2 == 0) {
      xaxis_text_colors = append(xaxis_text_colors, rep("red", length(sps)))  
    } else {
      xaxis_text_colors = append(xaxis_text_colors, rep("grey50", length(sps)))
    }
    color_counter = color_counter + 1
    
    # populate counts matrix with sp value (tumor%) or 0
    for (k in 1:length(sps)) {
      sp_data = sample_data[which(sample_data$SP==sps[k]),]
      sp_data = sp_data[grep("exonic", sp_data$region_type),]
      sp_genes = sort(unique(sp_data$region_name))
       
      for (l in 1:length(universe)) {
        if (universe[l] %in% sp_genes) { 
          #sample_heatmap[l,k] = sps[k] 
          sample_heatmap[l,k] = 1
        }
      }
    }
    
    # add sample heatmap to master heatmap
    master_heatmap = cbind(master_heatmap, sample_heatmap)
    
  } # end of sample

  # filter to non-singleton hits
  format_data = melt(master_heatmap)
  
  # plot mutation heatmap
  names(format_data) = c("variable", "sample", "value")
  base_size <- 12
  
  # adjust text size for patients with lots of mutated genes
  if (length(universe) > 80){
    resize= .2
  } else {
    resize= .6
  }
  
  yaxis_text_colors = universe %in% crc_genes
  yaxis_text_colors[yaxis_text_colors==FALSE] = "black"
  yaxis_text_colors[yaxis_text_colors==TRUE] = "red"
  
  (p <- ggplot(format_data, aes(x=sample, y=variable)) + geom_tile(aes(fill = value), colour="black") 
   + scale_fill_gradient(low = "white",high = "steelblue", limits=c(0,1)) 
   + theme_grey(base_size = base_size*resize) 
   + labs(x = "", y= "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
   + theme(legend.position = "none", axis.ticks = element_blank(), 
           axis.text.x = element_text(size=base_size*resize, angle=90, hjust=0, vjust=.5, colour=xaxis_text_colors),
           axis.text.y = element_text(colour=yaxis_text_colors)
           )
   + coord_fixed(ratio=1)
  )  
  ggsave(heatmap_img, p) 
  
} # end of sample groupings loop


