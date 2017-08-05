#!/usr/bin/python
#
# submit rsem transcriptome count
#
# usage: 
# python2 sumbitRsemExpression.py [fastq path] [sample name] [rsem ref] [proc]
# python2 submitRsemExpression.py \"../source_data/Sample_H_G$i/*.fastq\" Sample_H_G$i \"$HOME/references/rsem_resources/hg19\" 16"

import sys
import glob
import os
from subprocess import check_call, call
import re

# arguments
inputFastqFilePattern = sys.argv[1]
outputPrefix = sys.argv[2]
indexBaseDirectory = sys.argv[3]	#rsem reference $HOME/references/rsem_resources/hg19
procs = sys.argv[4]

# gets list of files based on the pattern passed at command line.
dirList = glob.glob(inputFastqFilePattern)
fastqFilesR1 = []
fastqFilesR2 = []
for fastq in dirList:
	if re.search("_R1_|_r1|\.1\.", fastq):
		fastqFilesR1.append(fastq)
	else:
		# assume single end read if doens't match for R1|r1
		fastqFilesR2.append(fastq)


# sort entries
fastqFilesR1 = sorted(fastqFilesR1)
fastqFilesR2 = sorted(fastqFilesR2)

print "fastq1 files."
print fastqFilesR1
print "fastq2 files."
print fastqFilesR2

# rsem-calculate-expression -p [procs] --paired-end [fastq1] [fastq2]
print "running rsem-calculate expression."
if len(fastqFilesR1) == 0:
	print "single end"
	cmd = 'rsem-calculate-expression -p ' + procs + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' ' + indexBaseDirectory + ' ' + outputPrefix
else:
	print "paired end"
	cmd = 'rsem-calculate-expression -p ' + procs + ' --paired-end ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + ' ' + indexBaseDirectory + ' ' + outputPrefix

print cmd
call(cmd, shell=True)	
