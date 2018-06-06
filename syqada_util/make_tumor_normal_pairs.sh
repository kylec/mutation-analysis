

PROJ=$1
cat pairs.txt | while read p t n; do echo -e "$t\t$n\t$t" >> $PROJ.tumor_normal.pairs; done
