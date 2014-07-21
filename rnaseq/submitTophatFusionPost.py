#!/usr/bin/python
import sys
import time
from popen2 import popen2
from sets import Set

runDirectory = sys.argv[1]
indexBaseDirectory = sys.argv[2]	# /RIS/home/fasaan/resources/tophat_resources/hg19/Homo_sapiens/UCSC/hg19
queue = sys.argv[3]
procs = sys.argv[4]
walltime = sys.argv[5]


# assemble cluster job
output, input = popen2('qsub')
input.write('#PBS -S /bin/bash\n')
input.write('#PBS -N ' + runDirectory + '\n')
input.write('#PBS -d ' + runDirectory +'\n')
input.write('#PBS -e ' + runDirectory + ' -o ' + runDirectory + '\n')
input.write('#PBS -q ' + queue + '\n')
input.write('#PBS -l procs=' + procs + ',walltime=' + walltime + '\n')

# commands for cluster job
input.write('tophat-fusion-post -p ' + procs + ' ')
input.write('--num-fusion-reads 2 ')
input.write('--num-fusion-pairs 2 ')
input.write('--num-fusion-both 8 ')
input.write('--output-dir tophatfusion_out_2_2_8 ')
input.write(indexBaseDirectory + '/Sequence/BowtieIndex/genome')

#finish script
input.close()
time.sleep(2)
	
#print('tophat -p ' + procs + ' ')
#print('--fusion-search ')
#print('-G ' + indexBaseDirectory + '/Annotation/Genes/genes.gtf ')
#print('-o ' + outputDirectory + ' ')
#print(indexBaseDirectory + '/Sequence/Bowtie2Index/genome ')
#print(','.join(fastqFilesR1) + '\n' + ','.join(fastqFilesR2))

