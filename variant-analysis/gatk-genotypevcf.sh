#!/usr/bin/env bash

#usage - sh gatk-genotypevcf.sh [output prefix] [build]

java8=/risapps/noarch/jdk/jdk1.8.0_45/bin/java

OUTPREFIX=$1
BUILD=$2

# get reference config
source ~/mutation-analysis/config.sh $BUILD

# ouputfiles
GVCFLIST=gvcf.list
RAWVCF=$OUTPREFIX.vcf
RECALSNPVCF=$OUTPREFIX.recalsnp.vcf
RECALSNPINDELVCF=$OUTPREFIX.recalsnpindel.vcf

# generate gvcf list
ls *.g.vcf > $GVCFLIST

# genotype gvcf files
if [ ! -f "$RAWVCF" ]; then 
echo "`date` Genotype gvcf"
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $REF --variant $GVCFLIST -o $RAWVCF -newQual
else
	echo "`date` $RAWVCF present"
fi

if [ "$?" != 0 ]; then echo "FAILURE: GenotypeGVCFs failed."; exit 1; fi

# build snp recal model
if [ ! -f "recalibrate_SNP.recal" ]; then
	echo "`date` snp recal model"
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF -input $RAWVCF -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE1_SNP -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R 
else 
	echo "`date` recalibrate_SNP.recal present."
fi

if [ "$?" != 0 ]; then echo "FAILURE: snp recal model failed."; exit 1; fi

# apply snp recal 
if [ ! -f "$RECALSNPVCF" ]; then
	echo "`date` apply snp recal."
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF -input $RAWVCF -mode SNP --ts_filter_level 99.5 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o $RECALSNPVCF
else 
	echo "`date` $RECALSNPVCF present"
fi

if [ "$?" != 0 ]; then echo "FAILURE: snp recal failed."; exit 1; fi

# build indel recal model
if [ ! -f "recalibrate_INDEL.recal" ]; then
	echo "`date` indel recal model"
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF -input $RECALSNPVCF -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL --maxGaussians 4 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -rscriptFile recalibrate_INDEL_plots.R 
else
	echo "`date` recalibrate_INDEL.recal present"
fi

if [ "$?" != 0 ]; then echo "FAILURE: indel recal model failed."; exit 1; fi

# apply indel recalibration
if [ ! -f "$RECALSNPINDELVCF" ]; then
	echo "`date` apply indel recal"
	$java8 -jar $HOME/bin/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF -input $RECALSNPVCF -mode INDEL --ts_filter_level 99.0 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -o $RECALSNPINDELVCF
else 
	echo "`date` $RECALSNPINDELVCF present"
fi

if [ "$?" != 0 ]; then echo "FAILURE: indel recal failed."; exit 1; fi

echo "`date` Done"
exit 0
