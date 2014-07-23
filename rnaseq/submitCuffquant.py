#!/usr/bin/python
import sys
import glob
import time
from popen2 import popen2
from sets import Set

inputFastqFilePattern = sys.argv[1]
runDirectory = sys.argv[2]
indexBaseDirectory = sys.argv[3]
queue = sys.argv[4]
nodes = sys.argv[5]
procs = sys.argv[6]
walltime = sys.argv[7]

sampleNames = set()
# gets list of files based on the pattern passed at command line.
dirList = glob.glob(inputFastqFilePattern)
for fname in dirList:
	# get name for output file
        filename = fname.split('/')[len(fname.split('/')) - 1]
	sampleName = fname.split('/')[len(fname.split('/')) - 2]
	print('Processing ' + sampleName)
	outputDirectory = sampleName + '_cqout'

	# assemble cluster job
	output, input = popen2('qsub')
	input.write('#PBS -S /bin/bash\n')
	input.write('#PBS -N ' + sampleName + '\n')
	input.write('#PBS -d ' + runDirectory +'\n')
	input.write('#PBS -e ' + runDirectory + ' -o ' + runDirectory + '\n')
	input.write('#PBS -q ' + queue + '\n')
	input.write('#PBS -l nodes=' + nodes + ':ppn=' + procs + ',walltime=' + walltime + '\n')

	# commands for cluster job
	input.write('cuffquant -p ' + procs + ' ')
	input.write('-o ' + outputDirectory + ' ')
	input.write('-b ' + indexBaseDirectory + '/Sequence/WholeGenomeFasta/genome.fa ')
	input.write(indexBaseDirectory + '/Annotation/Genes/genes.gtf ')
	input.write(fname + '\n')

	#finish script
	input.close()
	time.sleep(2)
	
	#print('tophat -p ' + procs + ' ')
	#print('--fusion-search ')
	#print('-G ' + indexBaseDirectory + '/Annotation/Genes/genes.gtf ')
	#print('-o ' + outputDirectory + ' ')
	#print(indexBaseDirectory + '/Sequence/Bowtie2Index/genome ')
	#print(','.join(fastqFilesR1) + '\n' + ','.join(fastqFilesR2))

