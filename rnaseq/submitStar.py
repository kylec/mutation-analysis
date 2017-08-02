#!/usr/bin/python
import sys
import glob
import os
from subprocess import check_call, call
import re

inputFastqFilePattern = sys.argv[1]
outputPrefix = sys.argv[2]
indexBaseDirectory = sys.argv[3] # ~/references/star_resources/	
procs = sys.argv[4]
gtf = '~/references/gencode.v19.annotation.nochr.gtf'

# 1st or 2nd pass switch
firstPass = sys.argv[5]

sampleNames = set()
# gets list of files based on the pattern passed at command line.
dirList = glob.glob(inputFastqFilePattern)
fastqFilesR1 = []
fastqFilesR2 = []
for fastq in dirList:
	if re.search("_R1_|_r1|\.1\.", fastq):
		fastqFilesR1.append(fastq)
	else:
		fastqFilesR2.append(fastq)

# sort entries
fastqFilesR1 = sorted(fastqFilesR1)
fastqFilesR2 = sorted(fastqFilesR2)

# add read groups
# ID:FCID_LANE PL:Illumina SM: Sample
readGroups = []
for i in fastqFilesR1:
	readGroups.append('ID:' + '_'.join(os.path.basename(i).split('_')[0:4]) + ' PL:Illumina SM:' + outputPrefix)
# if single-end reads
if len(fastqFilesR1) == 0:
	for i in fastqFilesR2:
		readGroups.append('ID:' + '_'.join(os.path.basename(i).split('_')[0:4]) + ' PL:Illumina SM:' + outputPrefix)

#print loaded fastq files
print "fastq1 files."
print fastqFilesR1
print "fastq2 files."
print fastqFilesR2

# STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN
print "running star."
if firstPass == '1':
	print 'Running first pass.'

	#cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --runThreadN ' + procs + '' + ' --outSAMtype None --outFileNamePrefix ' + outputPrefix
	cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --readFilesCommand zcat --runThreadN ' + procs + '' + ' --outSAMtype None --outFileNamePrefix ' + outputPrefix
else:
	print 'Running second pass.'
	#cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --runThreadN ' + procs + '' + ' --sjdbFileChrStartEnd combinedSJ.out.tab --limitSjdbInsertNsj 2000000 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFileNamePrefix ' + outputPrefix + ' --outSAMattrRGline ' + ' , '.join(readGroups)
	cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --readFilesCommand zcat --runThreadN ' + procs + '' + ' --sjdbFileChrStartEnd combinedSJ.out.tab --limitSjdbInsertNsj 2000000 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFileNamePrefix ' + outputPrefix + ' --outSAMattrRGline ' + ' , '.join(readGroups)
print cmd
call(cmd, shell=True)	
