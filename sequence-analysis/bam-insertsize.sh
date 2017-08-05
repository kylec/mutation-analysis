#get header
java8=/risapps/noarch/jdk/jdk1.8.0_45/bin/java
java=/risapps/noarch/jdk/jdk1.7.0_79/bin/java
PICARD_PATH=/risapps/noarch/picard/2.5.0/dist
sam=$1 
output=$2
TMPDIR=$PWD
MEM=16g

#sort sam
$java8 -Xmx10g -jar $PICARD_PATH/picard.jar CollectInsertSizeMetrics INPUT=$sam OUTPUT=$output.insertsize.txt H=$output.insertsize.pdf M=0.5

