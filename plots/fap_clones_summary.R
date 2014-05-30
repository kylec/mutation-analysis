# list of all mutated genes

library(ggplot2)
library(reshape2)

# read list of sample groupings
f = system("cat ~/Dropbox/lab_vilar/fap/samples_group.txt", intern=T)
tier_genes = read.table("~/Dropbox/lab_vilar/fap/expands/fap_mutect_tier_report", header=T, sep="\t")
tier_genes = tier_genes[c("tier", "chr", "hg19_pos")]

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

output_df = NULL

# get a list of genes for each samples
# create a whole set of genes
# each patient, each sample
for (i in 1:length(f)) {
  
  #### get universe genes by patient  
  s = gsub("\t", "*.tsv ", f[i])
  s = gsub("_[A-Z]", "", s)
  cmd = paste("awk -F \"\t\" \'$4==0\' ", s, "*.tsv | cut -f13-15 | grep exonic | egrep \"stop|nonsyn\" | cut -f2 | sort -u", sep="")
  mutated_genes = system(cmd, intern=T)
  universe = mutated_genes
  
  # get a list of samples (Vilar02_N) of a patient
  samples = strsplit(f[i], "\t")[[1]]
  
  # add genes of samples to a list
  for (j in 1:length(samples)) {
    orig_sample = strsplit(samples[j], "_")[[1]][1]
    fin = paste(orig_sample,"ann.tsv",sep=".")
    sample_data = read.table(fin, sep="\t", header=T, stringsAsFactors=F)
    
    head(sample_data)
    # each sample's sp, get genes
    sps = sort(unique(na.omit(sample_data$SP)), decreasing=TRUE)
    
    for (k in 1:length(sps)) {
      sample_data = merge(sample_data, tier_genes, by.x=c("chr", "startpos"), by.y=c("chr", "hg19_pos"))
      sp_data = sample_data[which(sample_data$SP==sps[k] & (sample_data$tier == "1"| sample_data$tier == "2"|
                                                              sample_data$tier == "3"|
                                                              sample_data$tier=="stop")),]
      genes = paste(sort(sp_data$genename), collapse=",")
      if (genes==""){
        genes="."
      }
      #genes = paste(sp_data$genename, sp_data$tier, sep="|", collapse=",")
      index = grep(orig_sample, fd$sample)
      if(length(index)>0){
        indel = as.character(fd[index,]$variable)
      } else {
        indel = "."
      }
      # print sample, sp, genes
      output = c(samples[j], length(sps), round(sps[k],2), genes, indel)
      cat(samples[j], length(sps), round(sps[k],2), genes, indel, "\n")
      output_df = rbind(output_df, output)
      
    }
    

  } # end of sample
  
  
} # end of sample groupings loop
colnames(output_df) = c("sample", "#clones", "Tumor%", "genes", "loh/indel")
write.table(output_df, "test.txt", quote=FALSE, row.names=FALSE, sep="\t")
