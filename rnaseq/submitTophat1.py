#!/usr/bin/python
import sys
import glob
import os
from subprocess import check_call, call
import re

inputFastqFilePattern = sys.argv[1]
outputDirectory = sys.argv[2]
indexBaseDirectory = sys.argv[3]	
procs = sys.argv[4]

# if theres sys.argv[8] run tophat mode not for fusionSearch 
# by default run tophat with fusionSearch
if len(sys.argv) > 4:
	fusionSearch = 0
else:
	fusionSearch = 1

sampleNames = set()
# gets list of files based on the pattern passed at command line.
dirList = glob.glob(inputFastqFilePattern)
print dirList
fastqFilesR1 = []
fastqFilesR2 = []
for fastq in dirList:
	if re.search('_R1_', fastq):
		fastqFilesR1.append(fastq)
	else:
		fastqFilesR2.append(fastq)

print fastqFilesR1
print "fastq2"
print fastqFilesR2

# tophat2 -o tophat_MCF7_1 -p 8 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM
	# commands for cluster job

if fusionSearch == 1:
	print "running with fusionSearch."
	cmd = 'tophat2 -p ' + procs + ' --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM -o ' + outputDirectory + ' ' + indexBaseDirectory + '/Sequence/BowtieIndex/genome ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2)
	
else:
	print "running without fusionSearch."
	cmd = 'tophat2 -p ' + procs + ' -G ' + indexBaseDirectory + '/Annotation/Genes/genes.gtf -o ' + outputDirectory + ' ' + indexBaseDirectory + '/Sequence/Bowtie2Index/genome ' + ','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2)
	
print cmd
call(cmd, shell=True)	
