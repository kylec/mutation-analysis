# Call mutect if .mutect output is not present
# Filter mutect output for PASS and COVERED positions. 
#
# Input: Tumor/normal BAMs
# Output: .mutect, .mutect.keep, .vcf, .pass.vcf
#
# usage
# 	sh mutect.sh [project name] [build]
# 	sh mutect.sh fap ucsc
#  
# kyle chang

# specifically use java 1.6 
JAVA=/scratch/rists/hpcapps/x86_64/jdk/jdk1.6.0_23/bin/java

PAT=$1
TSAM=$2
NSAM=$3
PROJ=$4
BUILD=$5
ANALYSISDIR=$6
BAMDIR=$ANALYSISDIR/$PROJ/bam/*.bam
SNVDIR=$ANALYSISDIR/$PROJ/mutect
COVDIR=$ANALYSISDIR/$PROJ/cov

# ucsc hg19
if [ "$BUILD" == "ucsc" ]; then
    REF=$ANALYSISDIR/references/ucsc.hg19.fasta
    COSMIC=$ANALYSISDIR/references/hg19_cosmic_v54_120711.vcf
    DBSNP=$ANALYSISDIR/references/dbsnp_137.hg19.vcf

# 1g ref or broad hg19 before ucsc hg19 release
# http://gatkforums.broadinstitute.org/discussion/2226/cosmic-and-dbsnp-files-for-mutect
elif [ "$BUILD" == "1g" ]; then
    REF=$ANALYSISDIR/references/Homo_sapiens_assembly19.fasta
    COSMIC=$ANALYSISDIR/references/b37_cosmic_v54_120711.vcf
    DBSNP=$ANALYSISDIR/references/dbsnp_137.b37.vcf
else 
    echo "ERROR: unknown build=$BUILD."
    exit 1
fi

# call mutect based on pair info
TBAM=`ls $BAMDIR | grep $TSAM`;
NBAM=`ls $BAMDIR | grep $NSAM`;
OUTPUTNAME=$PAT-$TSAM-$NSAM    
if [ -f "$SNVDIR/$OUTPUTNAME.mutect" ]; then
    echo "$SNVDIR/$OUTPUTNAME.mutect present."
else
    echo `date` mutect started.
    echo "java -Xmx4g -jar ~/bin/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $REF --dbsnp $DBSNP --cosmic $COSMIC --input_file:normal $NBAM --input_file:tumor $TBAM --out $SNVDIR/$OUTPUTNAME.mutect --coverage_file $COVDIR/$OUTPUTNAME.coverage.wig.txt --vcf $SNVDIR/$OUTPUTNAME.mutect.vcf"
    $JAVA -Xmx4g -jar ~/bin/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $REF --dbsnp $DBSNP --cosmic $COSMIC --input_file:normal $NBAM --input_file:tumor $TBAM --out $SNVDIR/$OUTPUTNAME.mutect --coverage_file $COVDIR/$OUTPUTNAME.coverage.wig.txt --vcf $SNVDIR/$OUTPUTNAME.mutect.vcf
    if [ "$?" == 0 ]; then
        # filter for pass mutations
        sed '1d' $SNVDIR/$OUTPUTNAME.mutect | awk -F"\t" '$35!="REJECT" && $10!="UNCOVERED"' > $SNVDIR/$OUTPUTNAME.mutect.keep
	grep -P "^#|PASS" $SNVDIR/$OUTPUTNAME.mutect.vcf > $SNVDIR/$OUTPUTNAME.mutect.pass.vcf
    else 
        echo "ERROR: mutect failed."
        exit 1
    fi
fi

echo `date` mutect done.
exit 0
