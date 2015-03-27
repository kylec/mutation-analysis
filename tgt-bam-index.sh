
# $1 = input bam
# $2 = target bed file
# $2 = output bam

intersectBed -abam $1 -b $2 > $3
samtools index $3 > $3.bai
