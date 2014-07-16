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

# scp files down
#for a in snv_only loh_sd_15 loh_mean_both_05 loh_mean_above_05; do 
for a in sampling-5 sampling-5-mean;do 
    #scp kchang3@D1prphaplotype0:/usr/local/epi/home/kchang3/fap/expands/$a/*sps* $a/
    echo $a
    cd $a
    
    for SAMPLE in `ls *.sps | cut -d. -f1`; do python /Users/kyle_air/Projects/mutation-analysis/mergeExpandsSNVAndVtoolsReport.py -i $SAMPLE.sps -o test > $SAMPLE.ann.tsv; done

    #### max sp per sample
    OUT="fap_expands_summary.txt"
    echo -e "Sample\tPurity\tClones\tSomatic\tLOH" > $OUT
    #for a in Vilar*.tsv; do 
    for a in *.sps; do 
        sample=`echo $a | cut -d. -f1`
        max_sp=`cut -f9 $a | sort -u | egrep -v "SP|NA" | sort -k1n | tail -1`
        num_sp=`cut -f9 $a | sort -u | egrep -v "SP|NA" | wc -l`
        somatic=`awk -F"\t" '$4==0' $a | wc -l`; loh=`awk -F"\t" '$4==1' $a | wc -l`
        echo -e "$sample\t$max_sp\t$num_sp\t$somatic\t$loh" >> $OUT
    done

    cd ..
done

# intersect pathway
cut -f1,4,5 /Users/kyle_air/Dropbox/lab_vilar/fap/patients.txt | grep -v Blood | while read sample type mut; do 
    echo $sample $type $mut
    a=`ls $sample.ann.tsv`
    pathway_list="$HOME/projects/fap/genes/crc-genes-pathway-nochr-vogel.bed"
    sed '1d' $a | awk '{FS=OFS="\t"; print $1, $2, $2, $0}' | cut -f1,2,3,6- > $a.bed
    intersectBed -a $a.bed -b $pathway_list -wo > $a.pathway
    
    sample=`echo $a | cut -d. -f1`
    max_sp=`cut -f9 $a | sort -u | egrep -v "SP|NA" | sort -k1n | tail -1`
    min_sp=`cut -f9 $a | sort -u | egrep -v "SP|NA" | sort -k1n | head -1`
    num_sp=`cut -f9 $a | sort -u | egrep -v "SP|NA" | wc -l`
    echo -e "Sample\t# Clones\t#Largest Clone\tSmallest Clone"
    echo -e "$sample\t$num_sp\t$max_sp\t$min_sp\n"
    echo -e "Clone\tGene\tMutationType\tCRC Pathway"
    grep exonic $a.pathway | cut -f10,15,16,22 | sort -rk1n 
    echo -e ""
done

## ad vogelstein to pathway
grep -w -f ../crc_genes/vogelstein.txt refseq_genes_061014.bed   | sort -u  > vogel.bed
sort -k1,1 -k2,2n vogel.bed > vogel.sorted.bed

cat crc-genes-pathway-nochr.bed vogel.bed > crc-genes-pathway-nochr-vogel.bed



#### filter evs
cut -f1,2 fap_mutect.report > all_var_coords
cut -f1,2 fap_mutect_novel_01_May23_181211.report > novel_01_var_coords
diff all_var_coords novel_01_var_coords  | grep "<" | cut -d" " -f2-  > excluded_var_coords
mkdir novel_01; cd novel_01
for a in `ls ../*combined`; do file=`echo $a | cut -d"/" -f2`;  grep -v -f excluded_var_coords $a > $file; done

# rename FAP files to vilar
for a in `ls FAP*`; do b=`echo $a | cut -d- -f2`; c=`echo $a | cut -d- -f2-`; mv $a $b-$c; done


###### varscan ######
for a in *dnacopy; do
    chr_count=`grep chr $a | wc -l`
    if [ "$chr_count" -gt "0" ]; then
        echo $a
        egrep -v "chrM|chrX|chrY" $a | sed 's/chr//g' > tmp && mv tmp $a
    fi
done

for a in *.called; do Rscript $HOME/mutation-analysis/dnacopy.R $a $a.dnacopy ; done
for a in *.copynumnber; do Rscript $HOME/mutation-analysis/dnacopy.R $a $a.dnacopy ; done

###### absolute ######
for a in `ls *pass.vcf | cut -d. -f1-3`; do 
    java -Xmx2g -jar ~/resources/snpEff/snpEff.jar -noStats -sequenceOntology -hgvs hg19 $a.vcf > $a.snpeff.vcf &
done

for a in `ls *pass.snpeff.vcf | cut -d. -f1`; do 
    tumor=`echo $a | cut -d- -f4-10`; normal=`echo $a | cut -d- -f11-17`; 
    perl ~/vcf2maf-master/vcf2maf.pl --input-snpeff $a.mutect.pass.snpeff.vcf --output-maf $a.maf --tumor-id $tumor --normal-id $normal
done

###### mutect ######
# filter mutect calls
for a in FAP*.snv; do echo $a; sed '1d' $a | awk -F"\t" '$35!="REJECT" && $10!="UNCOVERED"' > $a.keep; done

# filter mutect vcf calls
for a in `ls *.mutect.vcf | cut -d. -f1,2`; do grep -P "^#|PASS" $a.vcf > $a.pass.vcf; done

# filter mutect.keep by exome
capture_file="/usr/local/epi/home/kchang3/references/SeqCap_EZ_Exome_v3_primary.bed"
for a in `ls *.mutect.keep`; do
  echo $a; out=$a.exome;
  sed '1d' $a | awk '{FS=OFS="\t"; print $1, $2-1, $2, $0}' | cut -f1-3,6- > tmp
  head -n1 $a > $out
  intersectBed -a tmp -b $capture_file -wa | cut -f1,3- >> $out
  rm tmp
done

###### expands ######
cut_cmd="cut -d- -f1-3"
# make snv in put for (TCGA)
for a in *.keep ; do
    PAT=`echo $a | $cut_cmd`
    OUTFILE=$PAT.combined
    echo -e "chr\tstartpos\tAF_Tumor\tPN_B" > $OUTFILE
    sed '1d' $a | awk '{FS=OFS="\t"; print $1, $2, $22/($22+$21), 0}' >> $OUTFILE
done

# link inputs (TCGA)
mkdir expands; cd expands
ln -s ../mutect/*combined .
for a in ../varscan/*copynumber.dnacopy; do name=`basename $a | $cut_cmd`; echo $a $name.copynumber.dnacopy; done

# run expands
for a in *.combined; do sample=`echo $a | cut -d. -f1`; Rscript ~/mutation-analysis/expands.R $sample > $sample.log &  done