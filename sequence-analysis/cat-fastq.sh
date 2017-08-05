rm read1.fastq.gz; for a in `ls *R1*.gz`; do cat $a >> read1.fastq.gz; done
rm read2.fastq.gz; for a in `ls *R2*.gz`; do cat $a >> read2.fastq.gz; done
