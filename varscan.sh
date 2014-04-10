PROJ=$4
ANALYSISDIR=$HOME
BAMDIR=$ANALYSISDIR/$PROJ/bam/*.bam
OUTPUTDIR=$ANALYSISDIR/$PROJ/varscan
#REF=$ANALYSISDIR/references/ucsc.hg19.fasta
REF=$ANALYSISDIR/references/Homo_sapiens_assembly19.fasta

PAT=$1 
TSAM=$2
NSAM=$3

    #echo $PAT $TSAM $NSAM
	TBAM=`ls $BAMDIR | grep $TSAM`;
    NBAM=`ls $BAMDIR | grep $NSAM`;
    OUTPUTNAME=$OUTPUTDIR/$PAT-$TSAM-$NSAM
    echo $TBAM, $NBAM
    #echo `date` run samtools pileup.
	#samtools mpileup -q 1 -f $REF $NBAM $TBAM > $OUTPUTNAME.pileup

	# fix pileup fields
    if [ "$?" == 0 ]; then
        echo `date` fix pileup fileds
        #awk 'NF==9 && $4!=0' $OUTPUTNAME.pileup > $OUTPUTNAME.pileup.tmp && mv $OUTPUTNAME.pileup.tmp $OUTPUTNAME.pileup
    else 
        echo `date` pileup failed.
    fi

    if [ "$?" == 0 ]; then
        echo `date` varscan copynumber.
        java -jar ~/bin/VarScan.jar copynumber $OUTPUTNAME.pileup $OUTPUTNAME --mpileup 1
        #java -jar ~/bin/VarScan.jar copynumber $OUTPUTNAME.pileup $OUTPUTNAME --mpileup 1
	else
        echo `date` pileup failed.
    fi

    if [ "$?" == 0 ]; then
        echo `date` varscan copycaller.
        java -jar ~/bin/VarScan.jar copyCaller $OUTPUTNAME.copynumber --output-file $OUTPUTNAME.copynumber.called --output-homdel-file $OUTPUTNAME.copynumber.called.homdel
    else
        echo `date` copynumber failed.
    fi
