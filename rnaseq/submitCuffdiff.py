#!/usr/bin/python
import sys
from subprocess import call

mergedTranscripts = sys.argv[1]	# /RIS/home/fasaan/analysis/vilar_fap/rna_seq/clout/fap_hg19_cuffmerge/merged.gtf
conditions = sys.argv[2]
condition1List = sys.argv[3]
condition2List = sys.argv[4]
indexBaseDirectory = sys.argv[5]
outPrefix = sys.argv[6]
procs = sys.argv[7]

condition1Bams = []
condition1File = open(condition1List, 'r')
for line in condition1File:
	condition1Bams.append(line.strip())

condition2Bams = []
condition2File = open(condition2List, 'r')
for line in condition2File:
	condition2Bams.append(line.strip())


print("differential expression analysis for " + condition1List + " and " + condition2List + " (" + outPrefix + ")")

cmd = 'cuffdiff -u ' + mergedTranscripts + ' -b ' + indexBaseDirectory + '/Sequence/Chromosomes' + ' -p ' + procs + ' -o ' + outPrefix + ' -L ' + conditions + ' ' +  ','.join(condition1Bams) + ' ' + ','.join(condition2Bams)

print cmd

call(cmd, shell=True)

