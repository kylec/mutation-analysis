
RAWSAM=$1
SOBAM=$2

java8=/risapps/noarch/jdk/jdk1.8.0_45/bin/java
PICARD_PATH=/risapps/noarch/picard/2.5.0/dist

$java8 -Xmx10g -jar $PICARD_PATH/picard.jar SortSam INPUT=$RAWSAM OUTPUT=$SOBAM SORT_ORDER=coordinate
