# qsub wrapper to call mutect
#
# Input: project_name, build(ucsc, 1g), analysis_dir
#
# usage
# 	sh submitMutect.sh [project name] [build]
# 	sh submitMutect.sh fap ucsc
#  
# kyle chang


PROJ=$1
BUILD=$2
#analysis dir default $HOME or some other location
ANALYSISDIR=${3-$HOME} 
NODES=1
PROCS=8

if [ -f "$ANALYSISDIR/$PROJ/pairs.txt" ]; then
    # call mutect based on pair info
    cat $ANALYSISDIR/$PROJ/pairs.txt | while read PAT TSAM NSAM; do
        echo "submitted:  sh $SCRIPTSDIR/mutect.sh $PAT $TSAM $NSAM $PROJ $BUILD $ANALYSISDIR. qsub_args: $PAT $PAT.log $NODES $PROCS"
        q "sh $SCRIPTSDIR/mutect.sh $PAT $TSAM $NSAM $PROJ $BUILD $ANALYSISDIR" $PAT $PAT.log $NODES $PROCS
    done
else 
    echo "$ANALYSISDIR/$PROJ/pairs.txt missing."
    exit 1
fi
exit 0

