#!/usr/bin/python
import sys
import glob
import time
from popen2 import popen2
from sets import Set

mergedTranscripts = sys.argv[1]	# /RIS/home/fasaan/analysis/vilar_fap/rna_seq/clout/fap_hg19_cuffmerge/merged.gtf
conditions = sys.argv[2]
condition1List = sys.argv[3]
condition2List = sys.argv[4]
indexBaseDirectory = sys.argv[5]
outPrefix = sys.argv[6]
runDirectory = sys.argv[7]
queue = sys.argv[8]
nodes = sys.argv[9]
procs = sys.argv[10]
walltime = sys.argv[11]

condition1Bams = []
condition1File = open(condition1List, 'r')
for line in condition1File:
	condition1Bams.append(line.strip())

condition2Bams = []
condition2File = open(condition2List, 'r')
for line in condition2File:
	condition2Bams.append(line.strip())

# assemble cluster job
output, input = popen2('qsub')
input.write('#PBS -S /bin/bash\n')
input.write('#PBS -N ' + outPrefix + '\n')
input.write('#PBS -d ' + runDirectory +'\n')
input.write('#PBS -e ' + runDirectory + ' -o ' + runDirectory + '\n')
input.write('#PBS -q ' + queue + '\n')
input.write('#PBS -l nodes=' + nodes + ':ppn=' + procs + ',walltime=' + walltime + '\n')

print("differential expression analysis for " + condition1List + " and " + condition2List + " (" + outPrefix + ")")


# commands for cluster job
input.write('cuffdiff ')
input.write('-u ' + mergedTranscripts + ' ')
input.write('-b ' + indexBaseDirectory + '/Sequence/Chromosomes ')
input.write('-p ' + procs + ' ') 
input.write('-o ' + outPrefix + ' ')
input.write('-L ' + conditions + ' ') 
input.write(','.join(condition1Bams) + ' ')
input.write(','.join(condition2Bams) + '\n')

#finish script
input.close()

