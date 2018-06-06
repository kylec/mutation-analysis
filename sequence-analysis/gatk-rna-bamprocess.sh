#!/usr/bin/env bash

# usage
# sh gatk-rna-bamprocess [sample] [build(hg19/hg38)]

java8=/risapps/noarch/jdk/jdk1.8.0_45/bin/java
sam=$1 
BUILD=$2
TMPDIR=$PWD
MEM=16g

# source reference
source ~/mutation-analysis/config.sh $BUILD

# output
RAWBAM=$sam\Aligned.sortedByCoord.out.bam
HDBAM=$sam\Aligned.sortedByCoord.out.hd.bam
MDBAM=$sam\Aligned.sortedByCoord.out.hd.md.bam
STBAM=$sam\Aligned.sortedByCoord.out.hd.md.st.bam
RCBAM=$sam\Aligned.sortedByCoord.out.hd.md.st.rc.bam
RCTABLE=$sam\_recal_data.table
PRCTABLE=$sam\_post_recal_data.table
BQPLOT=$sam\_recalibration_plots.pdf


# run replace, markdup, split n trim 
if [ ! -f "$STBAM" ]; then
	echo "`date` Creating header."
	# create header
	header=$sam\header; samtools view -H $RAWBAM | egrep "@HD|@SQ|@PG" > $header; samtools view -H $RAWBAM | egrep "@RG" | uniq >> $header; samtools view -H $RAWBAM | egrep "@CO" >> $header;

	#replace header
	samtools reheader $sam\header $RAWBAM > $HDBAM

	echo "`date` mark dups."
	#mark duplicates
	$java8 -Xmx$MEM -jar $PICARDHOME/picard.jar MarkDuplicates I=$HDBAM O=$MDBAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$sam\_picard.metrics

	echo "`date` split and trim."
	#split n trim
	$java8 -Xmx$MEM  -jar ~/bin/GenomeAnalysisTK.jar -T SplitNCigarReads -R $REF -I $MDBAM -o $STBAM -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

	if [ "$?" != 0 ]; then echo "FAIL: mark dup, strim steps"; exit 1; fi
	
	# remove intermediates
	rm $HDBAM; rm $MDBAM
else 
	echo "`date` trim bam present."
fi

#### RUN REALIGN
#if [ ! -f "$RALIST" ]; then
#	echo "`date` create realign list"
#	# create target list for alignment (run once)
#	$java8 -Xmx$MEM -jar ~/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I $STBAM -known $MILLS -known $PHASE1_INDEL -o $RALIST
#	if [ "$?" != 0 ]; then echo "FAIL: create target list"; exit 1; fi
#else
#	echo "`date` $RALIST present."
#fi

#if [ ! -f "$RABAM" ]; then 
#	echo "`date` realign"
#	# realign
#	$java8 -Xmx$MEM -Djava.io.tmpdir=$TMPDIR -jar ~/bin/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $STBAM -targetIntervals $RALIST -known $MILLS -known $PHASE1_INDEL -o $RABAM
#	if [ "$?" != 0 ]; then echo "FAIL: realign"; exit 1; fi
#	# remove trim bam
#	rm $STBAM
#else
#	echo "`date` $RABAM present."
#fi

if [ ! -f "$RCBAM" ]; then
	echo "`date` base recal1."
	# base recal
	$java8 -Xmx$MEM -jar ~/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $STBAM -knownSites $DBSNP -knownSites $MILLS --knownSites $PHASE1_INDEL -o $RCTABLE
	if [ "$?" != 0 ]; then echo "FAIL: base recal1"; exit 1; fi

	echo "`date` base recal2."
	# base recal2 
	$java8 -Xmx$MEM -jar ~/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $STBAM -knownSites $DBSNP -knownSites $MILLS --knownSites $PHASE1_INDEL -BQSR $RCTABLE -o $PRCTABLE
	if [ "$?" != 0 ]; then echo "FAIL: base recal2"; exit 1; fi

	#echo "`date` base recal plot."
	# bq plots
	#$java -jar ~/bin/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF -before $RCTABLE -after $PRCTABLE -plots $BQPLOT
	#if [ "$?" != 0 ]; then echo "FAIL: base recal plot"; exit 1; fi

	echo "`date` base recal bam."
	# bq run
	$java8 -Xmx$MEM -jar ~/bin/GenomeAnalysisTK.jar -T PrintReads -R $REF -I $STBAM -BQSR $RCTABLE -o $RCBAM 
	if [ "$?" != 0 ]; then echo "FAIL: base recal bam"; exit 1; fi
else 
	echo "`date` $RCBAM present."
fi

