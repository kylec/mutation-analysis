#!/usr/bin/env bash

# usage
# sh gatk-haplotyper.sh [sample name] [bam] 

java8=/risapps/noarch/jdk/jdk1.8.0_45/bin/java
SAMPLE=$1
BAM=$2
DBSNP=$HOME/references/gatk-bundle/dbsnp_138.hg19.vcf
REF=$HOME/references/ucsc.hg19.fasta

$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF -I $BAM --emitRefConfidence GVCF --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 -o $SAMPLE.g.vcf --dbsnp $DBSNP
