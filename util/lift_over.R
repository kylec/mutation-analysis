library(rtracklayer)
#mutect_dna = readRDS("tcga_ucs_dnaseq_pass_mutect_snv.rds")
mutect_dna = data.frame(CHROM=c("chr1","chr1","chr1"), POS=c(1652516,100,1652516))
# liftover coordinates from hg38 to hg19
#import chain
chain = rtracklayer::import.chain("~/Downloads/hg38ToHg19.over.chain")
# create granges object for liftover
hg38_positions = with(mutect_dna, GRangesForUCSCGenome("hg38", CHROM, IRanges(mutect_dna$POS, mutect_dna$POS)))
# liftover
seqlevelsStyle(hg38_positions)
hg19_positions = liftOver(hg38_positions, chain)
length(hg38_positions) == length(hg19_positions)
# successfully lifted over positions
mutect_dna = mutect_dna[which(lengths(hg19_positions) > 0),]
length(which(lengths(hg19_positions) > 0)) == dim(mutect_dna)[1] 
# combine hg19 coordinates with mutect output
#mutect_dna = cbind(data.frame(CHROM=seqnames(unlist(hg19_positions)), POS=start(unlist(hg19_positions))), mutect_dna[,-(1:2)])
mutect_dna = cbind(data.frame(CHROM=seqnames(unlist(hg19_positions)), POS=start(unlist(hg19_positions))), mutect_dna)
mutect_dna$CHROM = as.character(mutect_dna$CHROM)

