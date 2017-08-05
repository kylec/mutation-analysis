#TODO: pass in tumor normal sample id from command line?
vcf=$1
#vcf="EBL11-EBL12.pass.vcf"
sample=`echo $vcf | cut -d. -f1`
polyp=`echo $vcf | cut -d- -f1`
normal=`grep -w $polyp ../pairs.txt | cut -f3`
echo $polyp $normal
prevep_vcf=$sample.prevep.vcf
vep_vcf=$sample.vep.vcf
grep "##" $vcf > $prevep_vcf
grep "CHROM" $vcf | cut -f1-10 >> $prevep_vcf
grep -v "#" $vcf | cut -f1-10 >> $prevep_vcf
#hla_file="../seq2hla/$polyp-ClassI.HLAgenotype4digits"
hla_file="../seq2hla/$normal-ClassI.HLAgenotype4digits"
program="NetMHCpan"
outdir="$polyp-result"
#annotate
if [ -f "$vep_vcf" ]; then
	echo "$vep_vcf present."
else
	perl ~/software/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl --input_file $prevep_vcf --format vcf --output_file $vep_vcf --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype --offline --species homo_sapiens_refseq --dir_plugins ~/.vep/Plugins/ --assembly GRCh37 --no_progress
fi

# get hla alleles

alleles=`cut -f2,4 $hla_file | grep -v Allele | sed 's/^/HLA-/' | sed "s/'//g" | tr '\t' "," | sed "s/,/,HLA-/g" | tr '\n' ',' | sed 's/,$//'`
# run pvacseq
#echo "pvacseq run $vep_vcf $polyp $alleles $program $outdir -e 8,9,10,11 --iedb-install-directory $HOME/software/IEDB/"
echo "pvacseq run $vep_vcf $polyp $alleles $program $outdir -e 8,9,10,11 --iedb-install-directory /usr/local/3rdparty/IEDB/"
pvacseq run $vep_vcf $polyp $alleles $program $outdir -e 8,9,10,11 --iedb-install-directory $HOME/software/IEDB/
