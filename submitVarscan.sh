PROJ=$1
BUILD=$2

cat $HOME/$PROJ/pairs.txt | while read PAT TSAM NSAM; do
    sh $HOME/scripts/varscan.sh $PAT $TSAM $NSAM $PROJ $BUILD &
done
