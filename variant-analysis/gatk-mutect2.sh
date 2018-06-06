# Call mutect if .mutect output is not present
# Filter mutect output for PASS and COVERED positions. 
#
# Input: Tumor/normal BAMs
# Output: .mutect, .mutect.keep, .vcf, .pass.vcf
#
# usage
# 	sh gatk-mutect2.sh [patient id] [tumor id] [normal id] [build] [analysis dir] [proj] [mutect dir] [bam dir]
#   sh gatk-mutect2.sh LS1 LS1_T LS1_N 1g ~/lynch_tumor rnaseq-hg19 mutect bam
#  
# kyle chang

JAVA=/risapps/noarch/jdk/jdk1.8.0_45/bin/java

PAT=$1
TSAM=$2
NSAM=$3
BUILD=$4
ANALYSISDIR=$5
PROJ=$6
MUTECTDIR=${7-"mutect"}
BAMDIRNAME=${8-"bam"}
BAMDIR=$ANALYSISDIR/$PROJ/$BAMDIRNAME/*.bam
SNVDIR=$ANALYSISDIR/$PROJ/$MUTECTDIR

# get reference config
source ~/mutation-analysis/config.sh $BUILD

# call mutect based on pair info
echo "$BAMDIR | grep -w $TSAM"
TBAM=`ls $BAMDIR | grep -w $TSAM`;
NBAM=`ls $BAMDIR | grep -w $NSAM`;
if [ ! -f "$TBAM" ] || [ ! -f "$NBAM" ]; then echo "FAILURE: t=$TBAM, n=$NBAM missing."; exit 1; fi
 
#OUTPUTNAME=$PAT-$TSAM-$NSAM    
OUTPUTNAME=$TSAM-$NSAM    
if [ -f "$SNVDIR/$OUTPUTNAME.vcf" ]; then
    echo "$SNVDIR/$OUTPUTNAME.vcf present."
else
    echo `date` mutect started.
    cmd="$JAVA -Djava.io.tmpdir=$SNVDIR -jar $HOME/bin/GenomeAnalysisTK.jar -T MuTect2 -R $REF --dbsnp $DBSNP --cosmic $COSMIC -I:tumor $TBAM -I:normal $NBAM -o $SNVDIR/$OUTPUTNAME.vcf"
		echo $cmd
		$cmd
    if [ "$?" != 0 ]; then echo "FAILURE: mutect failed."; exit 1; fi
fi

echo `date` mutect done.

if [ -f "$SNVDIR/$OUTPUTNAME.pass.vcf" ]; then
    echo "$SNVDIR/$OUTPUTNAME.pass.vcf present."
else
    echo `date` mutect pass started.
    sam=`echo $OUTPUTNAME | sed 's/-/\t/'`
    sed "s/TUMOR\tNORMAL/$sam/" $SNVDIR/$OUTPUTNAME.vcf | grep -P "^#|PASS"  > $SNVDIR/$OUTPUTNAME.pass.tmp
    mv $SNVDIR/$OUTPUTNAME.pass.tmp $SNVDIR/$OUTPUTNAME.pass.vcf
    if [ "$?" != 0 ]; then echo "FAILURE: mutect pass failed."; exit 1; fi
fi

echo `date` mutect done.

exit 0
