PROJ=$1
ANALYSISDIR=$HOME
BAMDIR=$ANALYSISDIR/$PROJ/bam/*.bam
SNVDIR=$ANALYSISDIR/$PROJ/mutect
COVDIR=$ANALYSISDIR/$PROJ/cov
#COSMIC=$ANALYSISDIR/references/chr_b37_cosmic_v54_120711.vcf
COSMIC=$ANALYSISDIR/references/b37_cosmic_v54_120711.vcf
#REF=$ANALYSISDIR/references/ucsc.hg19.fasta
REF=$ANALYSISDIR/references/Homo_sapiens_assembly19.fasta
#DBSNP=$ANALYSISDIR/references/dbsnp_132.hg19.vcf
DBSNP=$ANALYSISDIR/references/dbsnp_132_nochr.hg19.vcf

cat $ANALYSISDIR/$PROJ/pairs.txt | while read PAT TSAM NSAM; do
    #echo $PAT $TSAM $NSAM
	TBAM=`ls $BAMDIR | grep $TSAM`;
    NBAM=`ls $BAMDIR | grep $NSAM`;
    OUTPUTNAME=$PAT-$TSAM-$NSAM    
    echo $TBAM, $NBAM
	java -Xmx4g -jar ~/bin/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $REF --dbsnp $DBSNP --cosmic $COSMIC --input_file:normal $NBAM --input_file:tumor $TBAM --out $SNVDIR/$OUTPUTNAME.mutect --coverage_file $COVDIR/$OUTPUTNAME.coverage.wig.txt --vcf $SNVDIR/$OUTPUTNAME.mutect.vcf &
done
