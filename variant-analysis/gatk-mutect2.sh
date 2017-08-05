# Call mutect if .mutect output is not present
# Filter mutect output for PASS and COVERED positions. 
#
# Input: Tumor/normal BAMs
# Output: .mutect, .mutect.keep, .vcf, .pass.vcf
#
# usage
# 	sh gatk-mutect2.sh [patient id] [tumor id] [normal id] [build] [analysis dir] [proj] [bam dir] [mutect dir]
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

# ucsc hg19
if [ "$BUILD" == "ucsc" ]; then
    REF=$HOME/references/ucsc.hg19.fasta
    COSMIC=$HOME/references/hg19_cosmic_v54_120711.vcf
    DBSNP=$HOME/references/dbsnp_137.hg19.vcf

# 1g ref or broad hg19 before ucsc hg19 release
# http://gatkforums.broadinstitute.org/discussion/2226/cosmic-and-dbsnp-files-for-mutect
elif [ "$BUILD" == "1g" ]; then
    REF=$HOME/references/Homo_sapiens_assembly19.fasta
    COSMIC=$HOME/references/b37_cosmic_v54_120711.vcf
    DBSNP=$HOME/references/dbsnp_137.b37.vcf
elif [ "$BUILD" == "iot" ]; then
    REF=$HOME/references/hg19.fasta
    COSMIC=$HOME/references/hg19_iot_cosmic_v54_120711.vcf
    DBSNP=$HOME/references/dbsnp_137.iot.vcf
else 
    echo "ERROR: unknown build=$BUILD."
    exit 1
fi

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
    grep -P "^#|PASS" $SNVDIR/$OUTPUTNAME.vcf > $SNVDIR/$OUTPUTNAME.pass.tmp
    sam=`echo $OUTPUTNAME | sed 's/-/\t/'`
    sed "s/TUMOR\tNORMAL/$sam/" $SNVDIR/$OUTPUTNAME.vcf > $SNVDIR/$OUTPUTNAME.pass.tmp
    mv $SNVDIR/$OUTPUTNAME.pass.tmp $SNVDIR/$OUTPUTNAME.pass.vcf
    if [ "$?" != 0 ]; then echo "FAILURE: mutect pass failed."; exit 1; fi
fi

echo `date` mutect done.

exit 0
