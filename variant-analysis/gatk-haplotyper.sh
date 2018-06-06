#!/usr/bin/env bash

# usage
# sh gatk-haplotyper.sh [sample name] [build(hg38/hg19)]

SAMPLE=$1
BUILD=$2

BAM=`ls *recal.bam | grep -w $SAMPLE`

# get reference files
source ~/mutation-analysis/config.sh $BUILD

# run haplotyper
java -jar $HOME/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF -I $BAM --emitRefConfidence GVCF --genotyping_mode DISCOVERY -stand_call_conf 10 -o $SAMPLE.g.vcf --dbsnp $DBSNP -L $CAPTURE_FILE
