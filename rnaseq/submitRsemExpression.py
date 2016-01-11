#!/usr/bin/python
import sys
import glob
import os
from subprocess import check_call, call
import re

inputFastqFilePattern = sys.argv[1]
outputPrefix = sys.argv[2]
indexBaseDirectory = sys.argv[3]	#rsem reference $HOME/references/rsem_resources/hg19
procs = sys.argv[4]

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

#print "fastq1"
#print fastqFilesR1
#print "fastq2"
#print fastqFilesR2

# rsem-calculate-expression -p [procs] --paired-end [fastq1] [fastq2]
print "running rsem-calculate expression."
cmd = 'rsem-calculate-expression -p ' + procs + ' --paired-end ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' ' + indexBaseDirectory + ' ' + outputPrefix
	
print cmd
call(cmd, shell=True)	
