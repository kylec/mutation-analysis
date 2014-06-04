#### loh vs purity test on 3 samples
for sample in Vilar12 Vilar13 Vilar22; do 
    Rscript /Users/kyle_air/Projects/mutation-analysis/expands.R $sample > $sample.log
done


#### merge mutation and annotation
for SAMPLE in `ls *.sps | cut -d. -f1`; do python /Users/kyle_air/Projects/mutation-analysis/mergeExpandsSNVAndVtoolsReport.py -i $SAMPLE.sps -o test > $SAMPLE.ann.tsv; done
 
#### snv- only 
for a in `ls FAP*snv`; do b=`echo $a | cut -d- -f2`; mv $a $b.combined; done
for a in `ls FAP*dnacopy`; do b=`echo $a | cut -d- -f2`; mv $a $b.copynumber.dnacopy; done

#### count clones
dir=~/projects/fap/
cd $dir/expands
echo -e "sample\tdominant_clone\tnum_clones\tnonsyn\ton-target\traw"
for a in `ls *.ann.tsv | cut -d. -f1`; do 
    file=$a.ann.tsv; 
    mutect_exome=`ls $dir/mutect/$a*.exome`
    
    # highest tumor% clone; 
    sp=`cut -f9 $file | sort -u | sort -k1n | tail -1` 
    num_sp=`cut -f9 $file | sort -u | grep -v NA | wc -l`
    nonsyn=`grep $sp $file | grep exonic | cut -f13- | grep nonsyn | wc -l`
    ontarget=`grep -f<(cut -f1,2 $mutect_exome | sed 's/chr//g') $file | wc -l`
    raw=`sed '1d' $file | wc -l`
    echo -e "$a\t$sp\t$num_sp\t$nonsyn\t$ontarget\t$raw"
done


#### max sp per sample
for a in *.tsv; do sample=`echo $a | cut -d. -f1`; max_sp=`cut -f9 $a | sort -u | egrep -v "SP|NA" | sort -k1n | tail -1`; echo -e "$sample\t$max_sp"; done


#### filter evs
cut -f1,2 fap_mutect.report > all_var_coords
cut -f1,2 fap_mutect_novel_01_May23_181211.report > novel_01_var_coords
diff all_var_coords novel_01_var_coords  | grep "<" | cut -d" " -f2-  > excluded_var_coords
mkdir novel_01; cd novel_01
for a in `ls ../*combined`; do file=`echo $a | cut -d"/" -f2`;  grep -v -f excluded_var_coords $a > $file; done

# rename FAP files to vilar
for a in `ls FAP*`; do b=`echo $a | cut -d- -f2`; c=`echo $a | cut -d- -f2-`; mv $a $b-$c; done