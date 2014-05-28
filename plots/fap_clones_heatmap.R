# list of all mutated genes

library(ggplot2)
library(reshape2)

# read list of sample groupings
f = system("cat /Users/kyle_air/Dropbox/lab_vilar/fap/samples_group.txt", intern=T)

## add apc indel , loh data
fd=data.frame(variable = c(rep("APC-indel",2), rep("chr5-deletion", 5)),
              sample=c("Vilar15", "Vilar41", "Vilar12", "Vilar13", "Vilar16", "Vilar17", "Vilar56"),
              value=rep(1, 7))


# volgelstein gene list
vogelstein = system("cat ../crc_genes/vogelstein.txt", intern=T)
r_spondin = system("cat ../crc_genes/r_spondin.txt", intern=T)
hyper = system("cat ../crc_genes/tcga_hypermutated.txt", intern=T)
nonhyper = system("cat ../crc_genes/tcga_nonhypermutated.txt", intern=T)

#universe = mutated_genes
#heatmap_img = paste("fap_clones_mut_heatmap", "png", sep=".")

#patient_heatmap = NULL
#xaxis_text_colors = NULL
#color_counter = 0

# get a list of genes for each samples
# create a whole set of genes
# each patient, each sample
for (i in 1:length(f)) {
  
  #### get universe genes by patient
  patient_heatmap = NULL
  patient_indel_heatmap = NULL
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
    orig_sample = strsplit(samples[j], "_")[[1]][1]
    fin = paste(orig_sample,"ann.tsv",sep=".")
    sample_data = read.table(fin, sep="\t", header=T, stringsAsFactors=F)
    
    # each sample's sp, get genes
    sps = sort(unique(na.omit(sample_data$SP)), decreasing=TRUE)
    
    # create xlabel names by combing sample and sp 
    xaxis_text=c()
    for (q in 1:length(sps)) {
      xaxis_text = append(xaxis_text, paste(samples[j], round(sps[q],2), sep="_"))
    }
    
    # heat map matrix with tumor % 
    sample_heatmap <- matrix(0, nrow=length(universe), ncol=length(sps))
    colnames(sample_heatmap) <- xaxis_text
    rownames(sample_heatmap) = universe
    
    # alternate x axis label color to separate polyp groups
    if (color_counter %% 2 == 0) {
      xaxis_text_colors = append(xaxis_text_colors, rep("red", length(sps)))  
    } else {
      xaxis_text_colors = append(xaxis_text_colors, rep("grey50", length(sps)))
    }
    color_counter = color_counter + 1
    
    # populate counts matrix with sp value (tumor%) or 0
    # carry over mutated genes to lesser population
    carry_over_genes = NULL
    for (k in 1:length(sps)) {
      sp_data = sample_data[which(sample_data$SP==sps[k]),]
      sp_data = sp_data[grep("exonic", sp_data$region_type),]
      sp_genes = sort(unique(sp_data$region_name))
      sp_genes = append(sp_genes, carry_over_genes) 
      carry_over_genes = sp_genes
      
      for (l in 1:length(universe)) {
        if (universe[l] %in% sp_genes) { 
          sample_heatmap[l,k] = sps[k] 
          #sample_heatmap[l,k] = 1
        }
      }
    }
    
    # add sample heatmap to master heatmap
    patient_heatmap = cbind(patient_heatmap, sample_heatmap)
    
  } # end of sample

  # add indel/loh heatmap
  patient_indel_heatmap = matrix(0, nrow=3, ncol=length(colnames(patient_heatmap)))
  rownames(patient_indel_heatmap) = c("APC indel", "chr5 deletion","")
  colnames(patient_indel_heatmap) = colnames(patient_heatmap)
  for (y in 1:dim(fd)[1]) {
    index = grep(fd$sample[y], colnames(patient_heatmap))
    if (length(index) > 0) {
      patient_indel_heatmap[fd$variable[y], index] = 1
    }
  }
  
  patient_heatmap = rbind(patient_indel_heatmap, patient_heatmap)
  
  # annoates genes and genes list(vogelstein, hyper etc)
  yaxis_text = rownames(patient_heatmap)
  yaxis_text_colors = yaxis_text %in% c(vogelstein, r_spondin, hyper, nonhyper)
  yaxis_text_colors[yaxis_text_colors==FALSE] = "black"
  yaxis_text_colors[yaxis_text_colors==TRUE] = "red"
  for (y in 1:length(yaxis_text)) {
    if (yaxis_text[y] %in% vogelstein) {
      yaxis_text[y] = paste("*", yaxis_text[y])
    }
    if (yaxis_text[y] %in% r_spondin) {
      yaxis_text[y] = paste("^", yaxis_text[y])
    }
    if (yaxis_text[y] %in% hyper) {
      yaxis_text[y] = paste("o", yaxis_text[y])
    }
    if (yaxis_text[y] %in% nonhyper) {
      yaxis_text[y] = paste("x", yaxis_text[y])
    }
  }
  rownames(patient_heatmap) = yaxis_text
  genes_legend = data.frame(name=c("vogelstein", "r spondin", "hyper", "nonhyper"), value = c("*", "^", "o", "x"))
  
  # filter to non-singleton hits
  format_data = melt(patient_heatmap)
  
  # plot mutation heatmap
  names(format_data) = c("variable", "sample", "value")
  base_size <- 12
  
  # adjust text size for patients with lots of mutated genes
  if (length(universe) > 80){
    resize= .2
  } else {
    resize= .6
  } 
  
  (p <- ggplot(format_data, aes(x=sample, y=variable)) + geom_tile(aes(fill = value), colour="black") 
   + scale_fill_gradient(name="AF", low = "white",high = "steelblue", limits=c(0,1)) 
   + theme_grey(base_size = base_size*resize) 
   + labs(x = "", y= "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
   + theme(legend.position = "right", axis.ticks = element_blank(), 
           axis.text.x = element_text(size=base_size*resize, angle=90, hjust=0, vjust=.5, colour=xaxis_text_colors),
           axis.text.y = element_text(colour=yaxis_text_colors)
           )
   + coord_fixed(ratio=1)
  )  
  # todo: add geom_text for legend
  
  ggsave(heatmap_img, p) 
  
} # end of sample groupings loop

