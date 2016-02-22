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

# single end switch 
single_end = sys.argv[5]
# 1st or 2nd pass switch
firstPass = sys.argv[6]

sampleNames = set()
# gets list of files based on the pattern passed at command line.
dirList = glob.glob(inputFastqFilePattern)
#print dirList
fastqFilesR1 = []
fastqFilesR2 = []
for fastq in dirList:
	if re.search('_R1_', fastq):
		fastqFilesR1.append(fastq)
	else:
		fastqFilesR2.append(fastq)


# sort entries
fastqFilesR1 = sorted(fastqFilesR1)
fastqFilesR2 = sorted(fastqFilesR2)
#TODO: join fastqFiles as comma delimited strings for PE or SE to reduce if-else complexity during execution

# add read groups
# ID:FCID_LANE PL:Illumina SM: Sample
readGroups = []
for i in fastqFilesR1:
	readGroups.append('ID:' + '_'.join(os.path.basename(i).split('_')[0:4]) + ' PL:Illumina SM:' + outputPrefix)

#print "fastq1"
#print fastqFilesR1
#print "fastq2"
#print fastqFilesR2

# STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN
print "running star."
print single_end
if single_end == '1':
	cmd = ''
else:
	if firstPass == '1':
		print 'Running first pass.'
		# no gtf, bam output 
		#cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --readFilesCommand zcat --runThreadN ' + procs + '' + ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + outputPrefix + ' --outSAMattrRGline ' + ' , '.join(readGroups)
		# gtf, no bam output
		#cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --readFilesCommand zcat --runThreadN ' + procs + '' + ' --outSAMtype None --outFileNamePrefix ' + outputPrefix + ' --sjdbGTFfile ' + gtf
		# no gtf, no bam output
		cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --readFilesCommand zcat --runThreadN ' + procs + '' + ' --outSAMtype None --outFileNamePrefix ' + outputPrefix
	else:
		print 'Running second pass.'
		cmd = 'STAR --genomeDir ' + indexBaseDirectory + ' --readFilesIn ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' --readFilesCommand zcat --runThreadN ' + procs + '' + ' --sjdbFileChrStartEnd combinedSJ.out.tab --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFileNamePrefix ' + outputPrefix + ' --outSAMattrRGline ' + ' , '.join(readGroups)

print cmd
call(cmd, shell=True)	
