#### set project params
SCRIPT=$SCRIPTSDIR/rnaseq
QUEUE=long
NODES=1
PROCS=24:bigmem
WALLTIME=168:00:00
QUEUE=verylong
NODES=1
PROCS=24:bigmem
WALLTIME=1008:00:00
OUTPUTDIR=$HOME/fap/rnaseq-human
TOPHAT_REF=$HOME/references/tophat_resources/mmus10/Mus_musculus/UCSC/mm10
TOPHAT_REF=$HOME/references/tophat_resources/hg19/Homo_sapiens/UCSC/hg19
CUFFMERGE=fap_mm10_cuffmerge
CUFFMERGE=fap_hg19_cuffmerge

# (1) align rna with bowtie 1 for fusion search

# 1st batch of rnaseq
python $SCRIPT/submitTophat1.py \
"/RIS/home/scheet/projects/Vilar_FAP/sourcedata/set3-rnaseq-Hsap/**/*.fastq" \
"$OUTPUTDIR/thout" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME

# 2nd batch of rnaseq
python $SCRIPT/submitTophat1.py \
"/RIS/home/scheet/projects/Vilar_FAP/sourcedata/set8-rnaseq-Hsap/**/*.fastq.gz" \
"$OUTPUTDIR/thout" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME

# run tophat without fusionSearch, align directly to transcriptome
python $SCRIPT/submitTophat1.py \
"/RIS/home/scheet/projects/Vilar_FAP/sourcedata/set8-rnaseq-Hsap/**/*.fastq.gz" \
"$OUTPUTDIR/thout" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME nofusion

# (2) assemble expressed genes and transcripts
python $SCRIPT/submitCufflinks.py \
"$OUTPUTDIR/thout/*/accepted_hits.bam" \
"$OUTPUTDIR/clout" \
$QUEUE $NODES $PROCS $WALLTIME

# (3) create a file that lists transcripts for each sample
#     run this from the clout directory
ls **/transcripts.gtf > assemblies.txt

# (4) create a single merged transcriptome annotation for samples
python $SCRIPT/submitCuffmerge.py \
"$OUTPUTDIR/clout/assemblies.txt" \
$TOPHAT_REF \
$CUFFMERGE \
"$OUTPUTDIR/clout" \
$QUEUE $NODES $PROCS $WALLTIME

# (5) cuffquants
python $SCRIPT/submitCuffquant.py \
"$OUTPUTDIR/thout/*/accepted_hits.bam" \
"$OUTPUTDIR/cqout" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME

# (6) identify differentially expressed genes and transcripts
#for GROUP in large_intestine_polyp,large_instine_normal small_intestine_polyp,small_intestine_normal polyp,normal large_intestine,small_intestine 
FILE=bam
FILE=cxb
for GROUP in colon_polyp,colon_normal duodenum_polyp,duodenum_normal polyp,normal colon,duodenum; do
	GROUP1=`echo $GROUP | cut -d, -f1`
	GROUP2=`echo $GROUP | cut -d, -f2`
	
	python $SCRIPT/submitCuffdiff.py \
	"$OUTPUTDIR/clout/$CUFFMERGE/merged.gtf" \
	$GROUP \
	"$OUTPUTDIR/sample_groups/$GROUP1.$FILE.txt" \
	"$OUTPUTDIR/sample_groups/$GROUP2.$FILE.txt" \
	$TOPHAT_REF \
	$GROUP \
	"$OUTPUTDIR/cdout" \
	$QUEUE $NODES $PROCS $WALLTIME
done

# (7) identify fusions
python $SCRIPT/submitTophatFusionPost.py \
"$OUTPUTDIR/thfusion" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME
