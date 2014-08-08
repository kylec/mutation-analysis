Mutation Analysis
=================
## Project setup

```shell
PROJ=fap_ampliseq   
mkdir PROJ & cd $PROJ  
mkdir mutect bam varscan cov   
# Create samples.txt
# Link bams, vcf
cd bam; cut -f1,7 ../samples.txt  | while read sample bam; do echo $sample $bam; ln -s ../sourcedata/$bam.bai $sample.bam.bai; done
cut -f1,7 ../samples.txt  | while read sample bam; do echo $sample $bam; ln -s ../sourcedata/$bam $sample.bam; done
cut -f1,7 ../samples.txt  | while read sample bam; do echo $sample $bam; vcf=`echo $bam | sed 's/.bam/.vcf/'`; ln -s ../sourcedata/$vcf $sample.vcf; done
# Create pairs.txt - patient, tumor, normal

```

## Somatic mutation 
### Mutect
## Germline mutation
### GATK
## Copynumber
### Varscan2
## Clonality
### Chat
### Absolute
```shell
cat pairs.txt | while read pat tum nrm; do q "sh ~/mutation-analysis/absolute.sh $tum $tum $nrm fap ucsc $HOME" $tum $tum.log 1 2 1:00:00 short; done
grep -f samples.txt ../pairs.txt  | while read pat tum nrm ; do echo $pat $tum $nrm; q "sh ~/mutation-analysis/absolute.sh $tum $tum $nrm tcga_coadread ucsc /scratch/bcb/kchang3/Vilar_FAP" $tum $tum.log 1 2 1:00:00 short; done
```
### Expands
## Coverage
### Bedtools
```shell
for a in ../bam/*.bam; do echo $a; sample=`basename $a | cut -d. -f1`;  bedtools coverage -hist -abam $a -b /Users/kchang3/Analysis/references/IAD59895_182_Designed.bed > $sample.bam.hist.all; done
```

## RNAseq 
### Tuxedo Protocol
