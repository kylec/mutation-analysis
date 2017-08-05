#!/usr/bin/env bash
#
# usage
# sh samtools-callable-base.sh tumor.bam normal.bam tumor-normal.basecov.txt hg19

TBAM=$1
NBAM=$2
OUTPUTFILE=$3
BUILD=$4

if [ "$BUILD" == "hg19" ]; then
	echo "running with build $BUILD"
	REF=/rsrch2/ccp_rsch/kchang3/references/GRCh37-lite.fa
	samtools mpileup -q 1 -f /rsrch2/ccp_rsch/kchang3/references/GRCh37-lite.fa $TBAM $NBAM | awk -F "\t" '$4 > 9 && $7 > 9' | wc -l > $OUTPUTFILE
elif [ "$BUILD" == "hg38" ]; then
	echo "running with build $BUILD"
	REF=/rsrch2/epi/scheet/TEAM_ROOT/reference/genomics/Hsap/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
	LIST=/rsrch2/epi/scheet/TEAM_ROOT/reference/genomics/Hsap/GRCh38.d1.vd1/GRCh38.d1.vd1.chromosomes.bed
	samtools mpileup -q 1 -f $REF -l $LIST $TBAM $NBAM | awk -F "\t" '$4 > 9 && $7 > 9' | wc -l > $OUTPUTFILE
else
	echo "`date` wrong build=$BUILD"
	exit 1
fi

