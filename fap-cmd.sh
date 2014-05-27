

#### merge mutation and annotation
for SAMPLE in `ls *.sps | cut -d. -f1`; do python /Users/kyle_air/Projects/mutation-analysis/mergeExpandsSNVAndVtoolsReport.py -i $SAMPLE.sps -o test > $SAMPLE.ann.tsv; done
 
#### snv- only 
for a in `ls FAP*snv`; do b=`echo $a | cut -d- -f2`; mv $a $b.combined; done
for a in `ls FAP*dnacopy`; do b=`echo $a | cut -d- -f2`; mv $a $b.copynumber.dnacopy; done

#### count clones
cd ~/projects/fap/expands
echo -e "sample\tdominant_clone\tnum_clones\tnum_snv"
for a in `ls *.ann.tsv | cut -d. -f1`; do 
    file=$a.ann.tsv; 
    # highest tumor% clone; 
    sp=`cut -f9 $file | sort -u | sort -k1n | tail -1` 
    num_sp=`cut -f9 $file | sort -u | grep -v NA | wc -l`
    num_mut=`grep $sp $file | grep exonic | cut -f13- | grep nonsyn | wc -l`
    echo -e "$a\t$sp\t$num_sp\t$num_mut"
done


#### max sp per sample
for a in *.tsv; do sample=`echo $a | cut -d. -f1`; max_sp=`cut -f9 $a | sort -u | egrep -v "SP|NA" | sort -k1n | tail -1`; echo -e "$sample\t$max_sp"; done


#### filter out

Kyles-MacBook-Air:expands kyle_air$ cut -f1,2 fap_mutect.report > all_var_coords
Kyles-MacBook-Air:expands kyle_air$ cut -f1,2 fap_mutect_novel_01_May23_181211.report > novel_01_var_coords
diff all_var_coords novel_01_var_coords  | grep "<" | cut -d" " -f2-  > excluded_var_coords
mkdir novel_01; cd novel_01
for a in `ls ../*combined`; do file=`echo $a | cut -d"/" -f2`;  grep -v -f excluded_var_coords $a > $file; done