#### set project params

QUEUE=long
NODES=1
PROCS=24
WALLTIME=100:00:00
OUTPUTDIR=/RIS/home/scheet/projects/Vilar_FAP/rnaseq-human
OUTPUTDIR=/RIS/home/scheet/projects/Vilar_FAP/working/test/rna_seq
TOPHAT_REF=/RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19
CUFFMERGE=fap_hg19_cuffmerge

############ Human analysis ############

# (1) align rna with bowtie 1 for fusion search

# 1st batch of rnaseq
python $SCRIPTSDIR/submitTophat1.py \
"/RIS/home/scheet/projects/Vilar_FAP/sourcedata/set3-rnaseq-Hsap/**/*.fastq" \
"$OUTPUTDIR/thout" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME

# 2nd batch of rnaseq
python $SCRIPTSDIR/submitTophat1.py \
"/RIS/home/scheet/projects/Vilar_FAP/sourcedata/set8-rnaseq-Hsap/**/*.fastq.gz" \
"$OUTPUTDIR/thout" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME

# run tophat without fusionSearch, align directly to transcriptome
python $SCRIPTSDIR/submitTophat1.py \
"/RIS/home/scheet/projects/Vilar_FAP/sourcedata/set8-rnaseq-Hsap/**/*.fastq.gz" \
"$OUTPUTDIR/thout" \
$TOPHAT_REF \
$QUEUE $NODES $PROCS $WALLTIME nofusion

# (2) assemble expressed genes and transcripts
python $SCRIPTSDIR/submitCufflinks.py \
"$OUTPUTDIR/thout/*/accepted_hits.bam" \
"$OUTPUTDIR/clout" \
$QUEUE $NODES $PROCS $WALLTIME

# (3) create a file that lists transcripts for each sample
#     run this from the clout directory
ls **/transcripts.gtf > assemblies.txt

# (4) create a single merged transcriptome annotation for samples
python $SCRIPTSDIR/submitCuffmerge.py \
"$OUTPUTDIR/clout/assemblies.txt" \
$TOPHAT_REF \
$CUFFMERGE \
"$OUTPUTDIR/clout" \
$QUEUE $NODES $PROCS $WALLTIME

# (5) identify differentially expressed genes and transcripts
python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/clout/fap_hg19_cuffmerge/merged.gtf" \
"colon_polyp,colon_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/colon_polyp.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/colon_normal.txt" \
"/RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19" \
"colon_polyp_vs_colon_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/cdout" \
long 24 "200:00:00"

python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/clout/fap_hg19_cuffmerge/merged.gtf" \
"duodenum_polyp,duodenum_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/duodenum_polyp.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/duodenum_normal.txt" \
"/RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19" \
"duodenum_polyp_vs_duodenum_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/cdout" \
long 24 "200:00:00"

python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/clout/fap_hg19_cuffmerge/merged.gtf" \
"polyp,normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/polyp.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/normal.txt" \
"/RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19" \
"polyp_vs_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/cdout" \
long 24 "200:00:00"

python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/clout/fap_hg19_cuffmerge/merged.gtf" \
"colon,duodenum" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/colon.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/sample_groups/duodenum.txt" \
"/RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19" \
"colon_vs_duodenum" \ 
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/cdout" \
long 24 "200:00:00"

############ Mouse analysis ############

# (1) align rna-seq reads to the genome
python submitTophat.py \
"/RIS/home/scheet/projects/Vilar_FAP/sourcedata/set2-rnaseq-Mmus/**/*.fastq" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/thout" \
"/RIS/home/fasaan/resources/tophat_resources/mmus10/Mus_musculus/UCSC/mm10" \
long 24 "200:00:00"

# (2) assemble expressed genes and transcripts
python submitCufflinks.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/thout/**/accepted_hits.bam" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/clout" \
long 24 "200:00:00"

# (3) create a file that lists transcripts for each sample
#     run this from the clout directory
ls **/transcripts.gtf > assemblies.txt

# (4) create a single merged transcriptome annotation for samples
python submitCuffmerge.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/clout/assemblies.txt" \
"/RIS/home/fasaan/resources/tophat_resources/mmus10/Mus_musculus/UCSC/mm10" \
"fap_mm10_cuffmerge" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/clout" \
long 8 "24:00:00"

# (5) identify differentially expressed genes and transcripts
python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/clout/fap_mm10_cuffmerge/merged.gtf" \
"small_intestine_polyp,small_intestine_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/small_intestine_polyp.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/small_intestine_normal.txt" \
"/RIS/home/fasaan/resources/tophat_resources/mmus10/Mus_musculus/UCSC/mm10" \
"small_intestine_polyp_vs_small_intestine_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/cdout" \
long 24 "200:00:00"

python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/clout/fap_mm10_cuffmerge/merged.gtf" \
"large_intestine_polyp,large_intestine_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/large_intestine_polyp.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/large_intestine_normal.txt" \
"/RIS/home/fasaan/resources/tophat_resources/mmus10/Mus_musculus/UCSC/mm10" \
"large_intestine_polyp_vs_large_intestine_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/cdout" \
long 24 "200:00:00"

python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/clout/fap_mm10_cuffmerge/merged.gtf" \
"polyp,normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/polyp.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/normal.txt" \
"/RIS/home/fasaan/resources/tophat_resources/mmus10/Mus_musculus/UCSC/mm10" \
"polyp_vs_normal" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/cdout" \
long 24 "200:00:00"

python submitCuffdiff.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/clout/fap_mm10_cuffmerge/merged.gtf" \
"large_intestine,small_intestine" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/large_intestine.txt" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/sample_groups/small_intestine.txt" \
"/RIS/home/fasaan/resources/tophat_resources/mmus10/Mus_musculus/UCSC/mm10" \
"large_intestine_vs_small_intestine" \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq_mmus/cdout" \
long 24 "200:00:00"

######### Fusion ############
python submitTophatFusionPost.py \
"/RIS/home/fasaan/analysis/vilar_fap/rna_seq/thfusion" \
"/RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19" \
medium 24 "24:00:00"