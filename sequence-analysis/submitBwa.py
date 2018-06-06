#!/usr/bin/python

#example: python2.7 submitBwa.py dir/*fastq.gz sample_name ~/references/hg38.fasta 16
import sys
import glob
import os
from subprocess import check_call, call
import re

inputFastqFilePattern = sys.argv[1]
outputPrefix = sys.argv[2]
indexBaseDirectory = sys.argv[3] # ~/references/star_resources/	
procs = sys.argv[4]

print 'loading fastq...'
# gets list of files based on the pattern passed at command line.
dirList = glob.glob(inputFastqFilePattern)
print dirList
fastqFilesR1 = []
fastqFilesR2 = []
for fastq in dirList:
	if re.search("_R1_|_r1", fastq):
		fastqFilesR1.append(fastq)
	else:
		fastqFilesR2.append(fastq)


# sort entries
fastqFilesR1 = sorted(fastqFilesR1)
fastqFilesR2 = sorted(fastqFilesR2)
#TODO: join fastqFiles as comma delimited strings for PE or SE to reduce if-else complexity during execution

# add read groups
# @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1
# ID:FCID_LANE PL:Illumina SM: Sample
readGroup = '@RG\\tID:' + '_'.join(os.path.basename(fastqFilesR1[0]).split('_')[0:3]) + '\\tSM:' + outputPrefix + '\\tPL:illumina\\tLB:lib1\\tPU:unit1'

#print "fastq1"
print fastqFilesR1
#print "fastq2"
print fastqFilesR2
print readGroup
print 'running bwa...'

cmd = 'bwa mem -M -R \'' + readGroup + '\' -t 16 ' + indexBaseDirectory + ' ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' > ' + outputPrefix + '.sam'

print cmd
call(cmd, shell=True)	
