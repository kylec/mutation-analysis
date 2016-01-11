#!/bin/bash

# $1 = input bam name 
UNSORTED_BAM=$1.bam
SORTED_NAME=$1_namesort

echo "`date` sorting bam"
echo samtools sort -n $UNSORTED_BAM $SORTED_NAME
samtools sort -n $UNSORTED_BAM $SORTED_NAME
echo "`date` counting reads"
echo "htseq-count -f bam -s no -r name $SORTED_NAME.bam ~/references/tophat_resources/hg19/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > htseq-count.txt"
htseq-count -f bam -s no -r name $SORTED_NAME.bam ~/references/tophat_resources/hg19/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > htseq-count.txt
echo "`date` done"
