PROJ=$1
VCFPATH=$2
PAIRS=${3-"../pairs.txt"}

# usage
# ~/mutation-analysis/vtools_util/vtools_import.sh lynch "../mutect/*apcmmr.vcf" ../pairs.txt

if [ ! -f "$PROJ.proj" ]; then
	echo "`date` Init project $PROJ"
	vtools init $PROJ
else 
	echo "`date` $PROJ present"
fi

if [ ! -f "ANNOVAR.input" ]; then 
	echo "`date` load vcf"
	#for a in $VCFPATH; do vtools import --build hg19 --format control/mutect2_vcf.fmt $a; done
	cat $PAIRS | while read p t n; do 
		vcf=`ls $VCFPATH | grep $t-$n`
		echo $vcf
		if [ -f "$vcf" ]; then
			#vtools import --build hg19 --format ~/mutation-analysis/vtools_util/formats/mutect2_vcf.fmt $vcf --sample_name $p-$t-P $p-$n-N
			vtools import --build hg19 --format ~/mutation-analysis/vtools_util/formats/mutect2_vcf.fmt $vcf --sample_name $t $n
		else
			echo "WARNING: $vcf missing."
		fi 
	done
else
	echo "`date` annovar input present. assume vcf has been laoded. skip load vcf step."
fi

# vtools annotations
vtools use refGene-hg19_20130904
vtools use CancerGeneCensus-20130711 --linked_by refGene.name2
vtools use dbSNP-hg19_138
vtools use evs-6500 --as evs
vtools use thousandGenomes-hg19_20130502 --as thousandGenomes
vtools use ExAC-hg19_r0.2
vtools use dbNSFP
vtools use CosmicCodingMuts-v67_20131024 --as CosmicCodingMuts
vtools use CosmicMutantExport-v67_241013 --linked_by CosmicCodingMuts.COSMIC_ID
vtools use CosmicNonCodingVariants-v67_241013
vtools use ccdsGene-hg19_20130904
vtools use keggPathway-20110823 --linked_by ccdsGene.name
vtools export variant --output ANNOVAR.input --format ANNOVAR
ANNOVAR=/rsrch2/epi/scheet/TEAM_ROOT/3rdparty/no-arch/annovar-2016feb01/
perl $ANNOVAR/annotate_variation.pl -geneanno ANNOVAR.input -buildver hg19 $ANNOVAR/humandb/
vtools update variant --format ANNOVAR_exonic_variant_function --from_file ANNOVAR.input.exonic_variant_function --var_info mut_type function genename
vtools update variant --format ANNOVAR_variant_function --from_file ANNOVAR.input.variant_function --var_info region_type region_name

