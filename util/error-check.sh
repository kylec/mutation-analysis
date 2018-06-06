#error check

#gatk
for a in logs/*.log; do 
    finish=`tail -n1 $a | grep failing`
    if [ "$finish" != "" ]; then
        #echo $a
        tail -n1 $a
    else
        echo $a
    fi
done
# remove incomplete gatk file
find . -name "*.vcf" -size 0k | xargs rm
