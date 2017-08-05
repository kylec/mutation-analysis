#get header
java8=/risapps/noarch/jdk/jdk1.8.0_45/bin/java
java=/risapps/noarch/jdk/jdk1.7.0_79/bin/java
PICARD_PATH=/risapps/noarch/picard/2.5.0/dist
sam=$1 
TMPDIR=$PWD
MEM=16g
REF=$HOME/references/ucsc.hg19.fasta
#FIX<E change me back to .sam
RAWSAM=$sam.bam
SOBAM=$sam.sorted.bam
MDBAM=$sam.sorted.dedup.bam
RCBAM=$sam.sorted.dedup.recal.bam
RCBAI=$sam.sorted.dedup.recal.bai
HAPVCF=$sam.g.vcf

RALIST=$sam\_realignment_targets.list
RCTABLE=$sam\_recal_data.table
PRCTABLE=$sam\_post_recal_data.table
BQPLOT=$sam\_recalibration_plots.pdf
GOLDVCF1=$HOME/references/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf
GOLDVCF2=$HOME/references/gatk-bundle/1000G_phase1.indels.hg19.vcf
DBSNP=$HOME/references/gatk-bundle/dbsnp_138.hg19.vcf


# sort, dedup, index 
if [ ! -f "$MDBAM" ]; then
	echo "`date` Sort SAM."
	#sort sam
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar SortSam INPUT=$RAWSAM OUTPUT=$SOBAM SORT_ORDER=coordinate 

	echo "`date` mark dups."
	#mark duplicates
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar MarkDuplicates INPUT=$SOBAM OUTPUT=$MDBAM METRICS_FILE=$sam\_picard.metrics

	echo "`date` index."
	$java8 -Xmx10g -jar $PICARD_PATH/picard.jar BuildBamIndex INPUT=$MDBAM

	if [ "$?" != 0 ]; then echo "FAIL: sort,markdup,index"; exit 1; fi
	
	# remove intermediates
  #rm $RAWSAM	
else 
	echo "`date` dedup bam present."
fi

#exit 0

# recal
if [ ! -f "$RCBAM" ]; then
	echo "`date` base recal1."
	# base recal
	$java8 -Xmx$MEM -jar $HOME/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $MDBAM -knownSites $DBSNP -knownSites $GOLDVCF1 --knownSites $GOLDVCF2 -o $RCTABLE
	if [ "$?" != 0 ]; then echo "FAIL: base recal1"; exit 1; fi

	echo "`date` base recal2."
	# base recal2 
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $MDBAM -knownSites $DBSNP -knownSites $GOLDVCF1 --knownSites $GOLDVCF2 -BQSR $RCTABLE -o $PRCTABLE
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

if [ ! -f "$HAPVCF" ]; then
  
	if [ ! -f "$RCBAI" ]; then
		$java8 -Xmx10g -jar $PICARD_PATH/picard.jar BuildBamIndex INPUT=$RCBAM
	fi

	echo "`date` run haplotyper"
  $java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF -I $RCBAM --emitRefConfidence GVCF --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 -o $HAPVCF --dbsnp $DBSNP
  if [ "$?" != 0 ]; then echo "FAIL: haplotyper"; exit 1; fi
else
  echo "`date` $HAPVCF is present."
fi

echo "`date` Done"
exit 0
