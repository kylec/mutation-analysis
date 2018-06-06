#!/usr/bin/env bash

# gatk recal bam
# usage:
# sh gatk-dna-bamprocess.sh [SAMPLEple name] [build(hg19/hg38)]

java8=/risapps/noarch/jdk/jdk1.8.0_45/bin/java
SAMPLE=$1 
BUILD=$2
TMPDIR=$PWD
MEM=16g

# source reference
source ~/mutation-analysis/config.sh $BUILD

#FIX change me back to .SAMPLE
#RAWSAM=$SAMPLE.bam
RAWSAM=$SAMPLE.sam
SOBAM=$SAMPLE.sorted.bam
MDBAM=$SAMPLE.sorted.dedup.bam
RCBAM=$SAMPLE.sorted.dedup.recal.bam
RCBAI=$SAMPLE.sorted.dedup.recal.bai

RALIST=$SAMPLE\_realignment_targets.list
RCTABLE=$SAMPLE\_recal_data.table
PRCTABLE=$SAMPLE\_post_recal_data.table
BQPLOT=$SAMPLE\_recalibration_plots.pdf
ALIGN_SUM=$SAMPLE\_alignment_summary.txt


# sort, dedup, index 
if [ ! -f "$MDBAM" ]; then
	echo "`date` Sort SAM."
	#sort SAMPLE
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar SortSam INPUT=$RAWSAM OUTPUT=$SOBAM SORT_ORDER=coordinate 

	echo "`date` mark dups."
	#mark duplicates
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar MarkDuplicates INPUT=$SOBAM OUTPUT=$MDBAM METRICS_FILE=$SAMPLE\_picard.metrics

	echo "`date` index."
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar BuildBamIndex INPUT=$MDBAM

	if [ "$?" != 0 ]; then echo "FAIL: sort,markdup,index"; exit 1; fi
	
	# remove intermediates
  #rm $RAWSAM	
else 
	echo "`date` dedup bam present."
fi

# recal
if [ ! -f "$RCBAM" ]; then
	echo "`date` base recal1."
	# base recal
	$java8 -Xmx$MEM -jar $HOME/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $MDBAM -knownSites $DBSNP -knownSites $MILLS --knownSites $PHASE1_INDEL -o $RCTABLE
	if [ "$?" != 0 ]; then echo "FAIL: base recal1"; exit 1; fi

	echo "`date` base recal2."
	# base recal2 
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $MDBAM -knownSites $DBSNP -knownSites $MILLS --knownSites $PHASE1_INDEL -BQSR $RCTABLE -o $PRCTABLE
	if [ "$?" != 0 ]; then echo "FAIL: base recal2"; exit 1; fi

	echo "`date` base recal plot."
	# bq plots
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF -before $RCTABLE -after $PRCTABLE -plots $BQPLOT
	#if [ "$?" != 0 ]; then echo "FAIL: base recal plot"; exit 1; fi

	echo "`date` base recal bam."
	# bq run
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T PrintReads -R $REF -I $MDBAM -BQSR $RCTABLE -o $RCBAM 
	if [ "$?" != 0 ]; then echo "FAIL: base recal bam"; exit 1; fi
else 
	echo "`date` $RCBAM present."
fi

# alignment stats
if [ ! -f "$ALIGN_SUM" ]; then
	echo "`date` picard alignment summary."
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar CollectAlignmentSummaryMetrics R=$REF I=$RCBAM O=$ALIGN_SUM
	if [ "$?" != 0 ]; then echo "FAIL: alignment summary"; exit 1; fi
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar CollectInsertSizeMetrics I=$RCBAM O=$SAMPLE\_insertsize.txt H=$SAMPLE\_insertsize.pdf M=0.5
	if [ "$?" != 0 ]; then echo "FAIL: insert size metric"; exit 1; fi
else 
	echo "`date` $ALIGN_SUM present."
fi

echo "`date` Done"
exit 0
