PROJ=$1
BUILD=$2
ANALYSISDIR=$HOME
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

cat $ANALYSISDIR/$PROJ/pairs.txt | while read PAT TSAM NSAM; do
    #echo $PAT $TSAM $NSAM
	TBAM=`ls $BAMDIR | grep $TSAM`;
    NBAM=`ls $BAMDIR | grep $NSAM`;
    OUTPUTNAME=$PAT-$TSAM-$NSAM    
    echo $TBAM, $NBAM
	if [ -f "$SNVDIR/$OUTPUTNAME.mutect" ]; then
        echo "$SNVDIR/$OUTPUTNAME.mutect present."
    else 
        echo "$SNVDIR/$OUTPUTNAME.mutect not present."
        java -Xmx4g -jar ~/bin/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $REF --dbsnp $DBSNP --cosmic $COSMIC --input_file:normal $NBAM --input_file:tumor $TBAM --out $SNVDIR/$OUTPUTNAME.mutect --coverage_file $COVDIR/$OUTPUTNAME.coverage.wig.txt --vcf $SNVDIR/$OUTPUTNAME.mutect.vcf > $SNVDIR/$OUTPUTNAME.log &
    fi
done

exit 0
