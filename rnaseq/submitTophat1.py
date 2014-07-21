#!/usr/bin/python
import sys
import glob
import time
from popen2 import popen2
from sets import Set
import os

inputFastqFilePattern = sys.argv[1]
runDirectory = sys.argv[2]
indexBaseDirectory = sys.argv[3]	# /RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19
queue = sys.argv[4]
nodes = sys.argv[5]
procs = sys.argv[6]
walltime = sys.argv[7]

# if theres sys.argv[8] run tophat mode not for fusionSearch 
# by default run tophat with fusionSearch
if len(sys.argv) > 8:
	fusionSearch = 0
else:
	fusionSearch = 1

sampleNames = set()
# gets list of files based on the pattern passed at command line.
dirList = glob.glob(inputFastqFilePattern)
sampleToFilePrefix = {}
sampleToFastqMap = {}
for fname in dirList:
	# get name for output file
        filename = fname.split('/')[len(fname.split('/')) - 1]
	sampleName = fname.split('/')[len(fname.split('/')) - 2]
	sampleNames.add(sampleName)
	filenameTokens = filename.split('.')[0].split('_')
        fastqFilePrefix = filenameTokens[0] + '_' + filenameTokens[1]
	sampleToFilePrefix[sampleName] = fastqFilePrefix
	if sampleToFastqMap.get(sampleName) is None:
		sampleToFastqMap[sampleName] = []
	sampleToFastqMap[sampleName].append(fname)

# for each sample, generate the group level bam files
for sampleName in sampleNames:
	fastqFilePrefix = sampleToFilePrefix[sampleName]
	numFastqPairs = len(sampleToFastqMap[sampleName])/2
	outputDirectory = 'tophat_' + sampleName
	if os.path.exists(runDirectory + '/' + outputDirectory):
		print outputDirectory + 'exists.'
		continue
	print('Processing ' + str(numFastqPairs) + ' fastq pairs for ' + sampleName)
	fastqFilesR1 = []
	fastqFilesR2 = []
	for fastqFile in sampleToFastqMap[sampleName]:
		filename = fastqFile.split('/')[len(fname.split('/')) - 1].split('.')
		shortname = filename[0]
		#extension = filename[1]
		# correction for missing .gz extension in fastqFilesR2
		extension = '.'.join(filename[1:])
		shortnameTokens = shortname.split('_')
		directory = '/'.join(fastqFile.split('/')[:-1])
		if shortnameTokens[-2] == 'R2': continue		
		if shortnameTokens[-2] == 'R1':
			fastqFilesR1.append(fastqFile)
			shortnameTokens[-2] = 'R2'
			fastqFilesR2.append(directory + '/' + '_'.join(shortnameTokens) + '.' + extension)

	# assemble cluster job
	output, input = popen2('qsub')
	input.write('#PBS -S /bin/bash\n')
	input.write('#PBS -N ' + sampleName + '\n')
	input.write('#PBS -d ' + runDirectory +'\n')
	input.write('#PBS -e ' + runDirectory + ' -o ' + runDirectory + '\n')
	input.write('#PBS -q ' + queue + '\n')
	input.write('#PBS -l nodes=' + nodes + ':ppn=' + procs + ',walltime=' + walltime + '\n')

# tophat2 -o tophat_MCF7_1 -p 8 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM
	# commands for cluster job

	if fusionSearch == 1:
		print "running with fusionSearch."
		input.write('tophat2 -p ' + procs + ' ')
		input.write('--fusion-search ')
		input.write('--keep-fasta-order --bowtie1 --no-coverage-search -r 0 ')
        	input.write('--mate-std-dev 80 --fusion-min-dist 100000 --fusion-anchor-length 13 ')
       		input.write('--fusion-ignore-chromosomes chrM ')
		input.write('-o ' + outputDirectory + ' ')
		input.write(indexBaseDirectory + '/Sequence/BowtieIndex/genome ')
		input.write(','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + '\n')
	else:
		print "running without fusionSearch."
		# align speicific to gene transcript with -G 
        	input.write('tophat2 -p ' + procs + ' ')
        	input.write('--fusion-search ')
        	input.write('-G ' + indexBaseDirectory + '/Annotation/Genes/genes.gtf ')
        	input.write('-o ' + outputDirectory + ' ')
        	input.write(indexBaseDirectory + '/Sequence/Bowtie2Index/genome ')
        	input.write(','.join(fastqFilesR1) + ' ' + ','.join(fastqFilesR2) + '\n')
	
	#finish script
	input.close()
	time.sleep(2)
	
	#print('tophat -p ' + procs + ' ')
	#print('--fusion-search ')
	#print('-G ' + indexBaseDirectory + '/Annotation/Genes/genes.gtf ')
	#print('-o ' + outputDirectory + ' ')
	#print(indexBaseDirectory + '/Sequence/Bowtie2Index/genome ')
	#print(','.join(fastqFilesR1) + '\n' + ','.join(fastqFilesR2))

