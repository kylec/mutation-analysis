Mutation Analysis
=================
## Project setup

```shell
PROJ=fap_ampliseq   
mkdir PROJ & cd $PROJ  
mkdir mutect bam varscan cov   
# Create samples.txt
# Link bams, vcf
cd bam; cut -f1,7 ../samples.txt  | while read sample bam; do echo $sample $bam; ln -s ../sourcedata/$bam.bai $sample.bam.bai; done
cut -f1,7 ../samples.txt  | while read sample bam; do echo $sample $bam; ln -s ../sourcedata/$bam $sample.bam; done
cut -f1,7 ../samples.txt  | while read sample bam; do echo $sample $bam; vcf=`echo $bam | sed 's/.bam/.vcf/'`; ln -s ../sourcedata/$vcf $sample.vcf; done
# Create pairs.txt - patient, tumor, normal

```

## Somatic mutation 
### Mutect
### Varscan
```shell
# standard run in varscan dir
cat ../pairs.txt | while read pat tum nrm; do  
  q "sh varscan.sh $pat $tum $nrm fap ucsc $HOME" $pat $pat.log 1 16 168:00:00 long
done

#filter output by design
BEDFILE=~/references/IAD59895_182_Designed.bed
for VCF in *.snp.vcf; do
  OUTFILE=`basename $VCF | sed 's/.vcf/.tgt.vcf/'`
  bedtools intersect -a $VCF -b $BEDFILE > $OUTFILE
done

# generate p value segments for mergeSegments.pl
for a in `ls *.copynumber.called.dnacopy.seg | cut -d. -f1 | grep -v Vilar69`; do q "Rscript ~/mutation-analysis/dnacopy.R $a.copynumber.called $a.copynumber.called.dnacopy $a.copynumber.called.dnacopy_p" $a $a.log 1 8 10:00:00 medium; done
# format outputs after merge Segments
for sample in `ls *dnacopy_p | cut -d. -f1`; do 
  echo $sample
	input1=$sample.copynumber.called.dnacopy_p
	input2=$sample.copynumber.called.dnacopy
	merged=$sample.copynumber.called.dnacopy.merged
	awk '{FS=OFS="\t"; if ($2 ~ /#/ || $2 !~ /M|X|Y|GL|gl|hap/) { sub("chr","",$2); print $0}}' $input1 > tmp; mv tmp $input1
	awk '{FS=OFS="\t"; if ($1 ~ /#/ || $1 !~ /M|X|Y|GL|gl|hap/) { sub("chr","",$1); print $0}}' $input2 > tmp; mv tmp $input2
	perl ~/downloads/varscan-helper/mergeSegments.pl $input1 --ref-arm-sizes /scratch/bcb/kchang3/downloads/varscan-helper/cytoband.txt --output-basename $sample
	# format merged output 
	sed '1d' $sample.events.tsv | awk '{FS=OFS="\t"; print $1, $2, $3, $6, $4}' > $merged
done	
```
## Germline mutation
### GATK
## Copynumber
### Varscan2
## Clonality
### Chat
### Absolute
```shell
cat ../pairs.txt | while read pat tum nrm; do q "sh ~/mutation-analysis/absolute.sh $tum $tum $nrm fap ucsc $HOME varscan-adj absolute/dnacopy-adj" $tum $tum.log 1 2 1:00:00 short; done
grep -f samples.txt ../pairs.txt  | while read pat tum nrm ; do echo $pat $tum $nrm; q "sh ~/mutation-analysis/absolute.sh $tum $tum $nrm tcga_coadread ucsc /scratch/bcb/kchang3/Vilar_FAP" $tum $tum.log 1 2 1:00:00 short; done
```
### Expands
## Coverage
### Bedtools
```shell
for a in ../bam/*.bam; do echo $a; sample=`basename $a | cut -d. -f1`;  bedtools coverage -hist -abam $a -b /Users/kchang3/Analysis/references/IAD59895_182_Designed.bed > $sample.bam.hist.all; done
```

## RNAseq 
### Tuxedo Protocol
```shell
mkdir sourcedata sample_groups thout cqout clout cdout
#count fasta
a_path=/rsrch1/epi/scheet/PROJECTS/Vilar_Lynch/Project_EBF_LynchSyn_RNA48; b_path=.; for b in `ls $b_path`; do a_count=`ls $a_path/$b | wc -l`; b_count=`ls $b | wc -l`; echo -e "$b\t$a_count\t$b_count"; done
```
## Ampliseq
```shell
# torrent mutations
cat ../pairs.txt | while read pat tum nrm; do echo $tum $nrm; bedtools intersect -a $tum.vcf -b $nrm.vcf -v -header > $tum.somatic.vcf; done
for a in `ls *somatic.vcf` ; do echo $a;  sample=`echo $a | cut -d. -f1`; echo "python ~/mutation-analysis/adjustAf2vcf.py -i $a -b ../bam/$sample.bam > $sample.somatic.adjaf.vcf"; done
# import torrent mutations
cut -f1 ../pairs2.txt | while read tum; do echo $tum; vtools import --build hg19 --format control/torrent_fap_vcf.fmt ../torrent_somatic_low_stringency_15percent_adjaf/$tum.somatic.adjaf.vcf --sample_name $tum; done
less control/fap_torrent_report.header | vtools export variant --format control/fap_torrent_report.fmt --header - --output torrent.report --samples 'sample_name like "%EB%"'
python ~/mutation-analysis/addAf2Report.py -i torrent.report   > torrent.report.coding
awk -F"\t" '$71>0 || $1=="chr"' torrent.report.coding  > torrent.report.coding.afcutoff

# mutect mutations
BAMDIR=../../../0012-merge-and-filter/output
for PAT in `ls *.pass.vcf | cut -d. -f1`; do echo $PAT; python ~/mutation-analysis/adjustAf2vcf.py -i $PAT.mutect.pass.vcf -b $BAMDIR/$PAT.bam > $PAT.mutect.pass.adjaf.vcf
done
# import vtools
for a in `ls ~/fap_ampliseq_syqada/mutect-downsample-200/0001-mutect/output/*.adjaf.vcf`; do echo $a; tum=`basename $a | cut -d. -f1`; nrm=`grep $a ../pairs2.txt | cut -f3`; echo vtools import --build hg19 --format control/mutect_fap_vcf.fmt $a --sample_name $tum $nrm; done
less control/fap_torrent_report.header | vtools export variant --format control/fap_mutect_report.fmt --header - --output mutect.report --samples 'sample_name like "%-P"'

#vtools misc.
# replace sample name (not needed if you do vtools import by one sample, can be specified in vtools param)
cut -f1,2 ../samples.txt | while read curr prev; do echo $curr $prev; awk -v var1=$curr -v var2=$prev '{FS=OFS="\t"; if($1 ~ /#CHROM/){gsub(var2, var1, $10)}; print $0}' $curr.vcf > tmp; mv tmp $curr.vcf; done
# export novel 1% report
vtools select variant "(thousandGenomesEBI.EUR_AF_INFO is null or thousandGenomesEBI.EUR_AF_INFO < 0.01) and (evs.EuropeanAmericanMaf is null or 1.0*evs.EuropeanAmericanAltCount/(evs.EuropeanAmericanAltCount+evs.EuropeanAmericanRefCount) < 0.01)" -t variant_novel_01
less control/fap_torrent_report.header | vtools export variant_novel_01 --format control/fap_torrent_report.fmt --header - --output torrent_novel_01.report --samples 'sample_name like "%polyp"'
```

## Syqada
```shell
#syqada
# samples 
sed '1d'  samples.txt | cut -f2,3,7,8 | while read id pat bam char; do echo $pat-$id-$char; done > fap_ampliseq.samples
# pairs
cat pairs.txt  | while read pat tum nrm; do p=`grep $tum samples.txt | cut -f3`; t=`grep $tum samples.txt | awk '{OFS=""; print $3,"-",$2,"-",$8}'`; do t=`grep $tum samples.txt | awk '{OFS=""; print $3,"-",$2,"-",$8}'`; echo $p $t $n; done

#link bams
sed '1d'  samples.txt | cut -f2,3,7,8 | while read id pat bam char; do sample=$pat-$id-$char; echo $sample; ln -s ~/fap_ampliseq/sourcedata/$bam $sample.bam; done
```
