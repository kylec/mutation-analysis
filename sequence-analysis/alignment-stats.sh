# mapping statistics tophat
for a in `ls ../thout`; do p=`grep "concordant pair" ../thout/$a/align_summary.txt | cut -d" " -f1`; counts=`grep "Aligned pairs"  ../thout/$a/align_summary.txt | sed 's/Aligned pairs:  //g'`; echo -e "$a\t$counts\t$p"; done | sed 's/tophat_Sample_//g'
# mapping statistics star
grep "Uniquely mapped reads number" *Log.final.out | sed 's/Log.final.out://g' | sed 's/ *Uniquely mapped reads number *| *//g'
grep "Uniquely mapped reads %" *Log.final.out | sed 's/Log.final.out://g' | sed 's/ *Uniquely mapped reads % *| *//g'
