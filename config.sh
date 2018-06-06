BUILD=$1

#1. reference genome
if [ "$BUILD" == "hg19" ]; then
	REF=$HOME/references/ucsc.hg19.fasta
	COSMIC=$HOME/references/hg19_cosmic_v54_120711.vcf
	#DBSNP=$HOME/references/dbsnp_137.hg19.vcf
	DBSNP=$HOME/references/gatk-bundle/dbsnp_138.hg19.vcf
	MILLS=$HOME/references/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf
	PHASE1_INDEL=$HOME/references/gatk-bundle/1000G_phase1.indels.hg19.vcf
	PHASE1_SNP=$HOME/references/gatk-bundle/1000G_phase1.snps.high_confidence.hg19.vcf
	OMNI=$HOME/references/gatk-bundle/1000G_omni2.5.hg19.vcf
	HAPMAP=$HOME/references/gatk-bundle/hapmap_3.3.hg19.vcf
elif [ "$BUILD" == "grc37" ]; then
	# 1g ref or broad hg19 before ucsc hg19 release
	# http://gatkforums.broadinstitute.org/discussion/2226/cosmic-and-dbsnp-files-for-mutect
	REF=$HOME/references/Homo_sapiens_assembly19.fasta
	COSMIC=$HOME/references/b37_cosmic_v54_120711.vcf
	DBSNP=$HOME/references/dbsnp_137.b37.vcf
	MILLS=$HOME/references/gatk-bundle/Mills_and_1000G_gold_standard.indels.b37.vcf
	PHASE1_INDEL=$HOME/references/gatk-bundle/1000G_phase1.indels.b37.vcf
	PHASE1_SNP=~/references/gatk-bundle/1000G_phase1.snps.high_confidence.b37.vcf
	OMNI=$HOME/references/gatk-bundle/1000G_omni2.5.b37.vcf
	HAPMAP=$HOME/references/gatk-bundle/hapmap_3.3.b37.vcf
elif [ "$BUILD" == "iot" ]; then
	REF=$HOME/references/hg19.fasta
	COSMIC=$HOME/references/hg19_iot_cosmic_v54_120711.vcf
	DBSNP=$HOME/references/dbsnp_137.iot.vcf
elif [ "$BUILD" == "hg38" ]; then
	REF=$HOME/references/hg38/gatk-bundle/Homo_sapiens_assembly38.fasta
	#COSMIC=
	DBSNP=$HOME/references/hg38/gatk-bundle/dbsnp_146.hg38.vcf.gz
	MILLS=$HOME/references/hg38/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
	PHASE1_INDEL=$HOME/references/hg38/gatk-bundle/Homo_sapiens_assembly38.known_indels.vcf.gz
	PHASE1_SNP=$HOME/references/hg38/gatk-bundle/1000G_phase1.snps.high_confidence.hg38.vcf
	OMNI=$HOME/references/hg38/gatk-bundle/1000G_omni2.5.hg38.vcf
	HAPMAP=$HOME/references/gatk-bundle/hapmap_3.3.hg38.vcf
else
    echo "ERROR: unknown build=$BUILD."
    exit 1
fi

#2. capture file
CAPTURE_FILE=$HOME/references/hg38/S30409818_Padded_hg38.bed
#CAPTURE_FILE=$HOME/references/hg38/test-capture.bed

