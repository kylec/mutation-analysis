# Annotatae mutect vcf with snpEff
# Convert to MAF
# run absolute
#
# usuage
#	sh absolute.sh patient tumor normal ucsc $HOME
#
# Kyle Chang

PAT=$1
TSAM=$2
NSAM=$3
PROJ=$4
BUILD=$5
ANALYSISDIR=$6


SNVDIR=$ANALYSISDIR/$PROJ/mutect

if [ "$BUILD" == "ucsc" ]; then
    REF="hg19"
elif [ "$BUILD" == "1g" ]; then
    REF="GRCh37.75"
else
    echo "ERROR:build==$BUILD not recognized."
    exit 1
fi

# check mutect and cnv files exist
#VCF=`ls $SNVDIR/*.pass.exome.vcf | grep $PAT`
#VCF=`ls $SNVDIR/*.pass.vcf | grep -w $TSAM`
VCF=`ls $SNVDIR/*.vcf | grep -v pass | grep -w $TSAM`
if [ "$VCF" == "" ]; then
    echo "ERROR: missing file for $PAT, vcf=$VCF, cnv=$CNV."; exit 1
fi

#VCF_EFF=`echo $VCF | sed 's/pass/pass.snpeff/'`
VCF_EFF=`echo $VCF | sed 's/vcf/snpeff.vcf/'`
#MAF=`echo $VCF | sed 's/pass.vcf/pass.snpeff.maf/'`
MAF=`echo $VCF | sed 's/vcf/snpeff.maf/'`

# run SNPEFF
if [ -f "$VCF_EFF" ]; then
    echo "`date` $VCF_EFF present."
else 
    echo "`date` Run snpeff on vcf"
    echo "java -Xmx4G -jar $SNPEFF_HOME/snpEff.jar -noStats -sequenceOntology -hgvs $REF $VCF > $VCF_EFF"
    java -Xmx8G -jar $SNPEFF_HOME/snpEff.jar -noStats -sequenceOntology -hgvs $REF $VCF > $VCF_EFF
fi

# VCF to MAF
if [ "$?" == 0 ]; then
    if [ -f "$MAF" ]; then		
        echo "`date` $MAF present."
    else
        # remove chrX,Y,M, gl
        echo "`date` Remove chrX,Y,M, gl and chr prefix"
        awk '{FS=OFS="\t"; if ($1 ~ /#/ || $1 !~ /M|X|Y|GL|gl|hap/) { print $0}}' $VCF_EFF | sed 's/^chr//g' > $VCF_EFF.tmp && mv $VCF_EFF.tmp $VCF_EFF
        echo "`date` Convert vcf to maf."
        echo "perl $VCF2MAF_HOME/vcf2maf.pl --input-snpeff $VCF_EFF --output-maf $MAF --tumor-id $TSAM --normal-id $NSAM"
        perl $VCF2MAF_HOME/vcf2maf.pl --input-snpeff $VCF_EFF --output-maf $MAF --tumor-id $TSAM --normal-id $NSAM
    fi
else
    echo "ERROR: snpeff error."; exit 1
fi

if [ "$?" == 0 ]; then
    echo "`date` Done"
else 
    echo "ERROR: VCF2MAF failed."; exit 1 
fi
