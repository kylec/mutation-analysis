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

## Ampliseq
```shell
# mutations
cat ../pairs.txt | while read pat tum nrm; do echo $tum $nrm; bedtools intersect -a $tum.vcf -b $nrm.vcf -v -header > $tum.somatic.vcf; done
# replace sample name (not needed if you do vtools import by one sample, can be specified in vtools param)
cut -f1,2 ../samples.txt | while read curr prev; do echo $curr $prev; awk -v var1=$curr -v var2=$prev '{FS=OFS="\t"; if($1 ~ /#CHROM/){gsub(var2, var1, $10)}; print $0}' $curr.vcf > tmp; mv tmp $curr.vcf; done
#  import
cut -f1,2,3,5 ../samples.txt | while read sample id pat type; do vcf=`ls ../torrent/*.somatic.vcf | grep $sample`; if [ "$vcf" != "" ]; then vtools import --build hg19 --format control/torrent_fap_vcf.fmt $vcf --sample_name $pat-$id-$type; fi; done
# export vtools report
less control/fap_torrent_report.header | vtools export variant --format control/fap_torrent_report.fmt --header - --output torrent.report --samples 'sample_name like "%polyp"'
# export novel 1% report
vtools select variant "(thousandGenomesEBI.EUR_AF_INFO is null or thousandGenomesEBI.EUR_AF_INFO < 0.01) and (evs.EuropeanAmericanMaf is null or 1.0*evs.EuropeanAmericanAltCount/(evs.EuropeanAmericanAltCount+evs.EuropeanAmericanRefCount) < 0.01)" -t variant_novel_01
less control/fap_torrent_report.header | vtools export variant_novel_01 --format control/fap_torrent_report.fmt --header - --output torrent_novel_01.report --samples 'sample_name like "%polyp"'
```
## Syqada
```shell
#syqada
# samples 
sed '1d'  samples.txt | cut -f2,3,7,8 | while read id pat bam char; do echo $pat-$id-$char; done > fap_ampliseq.samples
# pairs
cat pairs.txt  | while read pat tum nrm; do p=`grep $tum samples.txt | cut -f3`; t=`grep $tum samples.txt | awk '{OFS=""; print $3,"-",$2,"-",$8}'`; do t=`grep $tum samples.txt | awk '{OFS=""; print $3,"-",$2,"-",$8}'`; echo $p $t $n; done

#link bams
sed '1d'  samples.txt | cut -f2,3,7,8 | while read id pat bam char; do sample=$pat-$id-$char; echo $sample; ln -s ~/fap_ampliseq/sourcedata/$bam $sample.bam; done
```
