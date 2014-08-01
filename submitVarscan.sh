## qsub wrapper for submitting varscan
## 
## Input: project_name, build(ucsc,1g), analysis_dir
##
## Usage:
## 	sh submitVarscan.sh fap ucsc $HOME
##
## kyle chang 

PROJ=$1
BUILD=$2
#analysis dir default home or another location
ANALYSISDIR=${3-$HOME}
NODES=1
PROCS=8

if [ -f "$ANALYSISDIR/$PROJ/pairs.txt" ]; then
    # call mutect based on pair info
    cat $ANALYSISDIR/$PROJ/pairs.txt | while read PAT TSAM NSAM; do
        OUTPUTFILE=`ls $ANALYSISDIR/$PROJ/varscan/*.copynumber | grep $TSAM`
        if [ "$OUTPUTFILE" == "" ]; then
            echo "submitted:  sh $SCRIPTSDIR/varscan.sh $PAT $TSAM $NSAM $PROJ $BUILD $ANALYSISDIR. qsub_args: $PAT $PAT.log $NODES $PROCS"
            q "sh $SCRIPTSDIR/varscan.sh $PAT $TSAM $NSAM $PROJ $BUILD $ANALYSISDIR" $PAT $PAT.log $NODES $PROCS
        else
            echo "$OUTPUTFILE exists."
        fi
    done
else
    echo "$ANALYSISDIR/$PROJ/pairs.txt missing."
    exit 1
fi
exit 0

