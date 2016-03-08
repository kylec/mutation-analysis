#############################
# usage
# sh script.sh [project_dir] [bam dir] [snp list]
# e.g. sh script.sh "~/Analysis/fap/rnaseq-hg19" "~/Analysis/fap/rnaseq-hg19/*.rc.bam" input.txt
##############################

# input arguments
PROJDIR=$1
BAMDIR=$2
INPUT=$3

# print output
OUTPUT=cmd.txt
echo snapshotDirectory $PROJDIR/igv > $OUTPUT
echo genome hg19 >> $OUTPUT

# read input
sed '1d' $PROJDIR/igv/$INPUT | cut -f1-3 | while read SAMPLE CHR START; do
#cat $PROJDIR/igv/$INPUT | cut -f1-3 | while read SAMPLE CHR START; do
	CHR=`echo $CHR | sed 's/^chr//g'`	
	# assume $SAMPLE is a polyp, use it to find patient id and normal sample id
    PATIENT=`grep -w $SAMPLE $PROJDIR/pairs.txt | cut -f1`
    TUMOR=`grep -w $SAMPLE $PROJDIR/pairs.txt | cut -f2`
    NORMAL=`grep -w $SAMPLE $PROJDIR/pairs.txt | cut -f3`
	
	# zoom out
	RANGE1=$((START-150))
	RANGE2=$((START+150))

	#echo $SAMPLE
	#echo "grep -w $SAMPLE $PROJDIR/samples.txt | cut -f3"
	#echo "$PATIENT, $TUMOR, $NORMAL"
	echo "ls $BAMDIR | grep -w $TUMOR"
	TBAM=`ls $BAMDIR | grep -w $TUMOR`
	NBAM=`ls $BAMDIR | grep -w $NORMAL`

	echo "tbam=$TBAM nbam=$NBAM"
	
	if [ -e "$PATIENT-$SAMPLE-polyp-$CHR-$START.jpg" ]; then
		echo exists.
	else 
		echo $PATIENT $TUMOR $NORMAL $START $RANGE1 $RANGE2
		# print igv command
		echo new >> $OUTPUT
		echo load $TBAM >> $OUTPUT
		echo goto chr$CHR:$START-$START >> $OUTPUT
		#echo snapshot $PATIENT-$SAMPLE-polyp-$CHR-$START.jpg >> $OUTPUT
 		echo collapse $TBAM >> $OUTPUT
		echo snapshot $PATIENT-$SAMPLE-polyp-$CHR-$START-collapse.jpg >> $OUTPUT
		echo goto chr$CHR:$RANGE1-$RANGE2 >> $OUTPUT
   	 	echo snapshot $PATIENT-$SAMPLE-polyp-$CHR-$START-zoomout-collapse.jpg >> $OUTPUT
		#echo expand $TUMOR.bam >> $OUTPUT
		#echo snapshot $PATIENT-$SAMPLE-polyp-$CHR-$START-zoomout.jpg >> $OUTPUT
		
		echo new >> $OUTPUT
		echo load $NBAM >> $OUTPUT
		echo goto chr$CHR:$START-$START >> $OUTPUT
		#echo snapshot $PATIENT-$SAMPLE-normal-$CHR-$START.jpg >> $OUTPUT
		echo collapse $NBAM >> $OUTPUT
		echo snapshot $PATIENT-$SAMPLE-normal-$CHR-$START-collapse.jpg >> $OUTPUT
		echo goto chr$CHR:$RANGE1-$RANGE2 >> $OUTPUT
		echo snapshot $PATIENT-$SAMPLE-normal-$CHR-$START-zoomout-collapse.jpg >> $OUTPUT
		#echo expand $NORMAL.bam >> $OUTPUT
    	#echo snapshot $PATIENT-$SAMPLE-normal-$CHR-$START-zoomout.jpg >> $OUTPUT
	fi	
done
#echo exit >> $OUTPUT

# execute batch file
#echo execute batch file...
#java -Xmx4000m \
#    -Dapple.laf.useScreenMenuBar=true \
#    -Djava.net.preferIPv4Stack=true \
#    -jar /Users/kchang3/Downloads/resources/igv/igv.jar -b $OUTPUT
