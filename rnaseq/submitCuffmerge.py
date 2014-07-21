#!/usr/bin/python
import sys
import glob
import time
from popen2 import popen2
from sets import Set

assembliesFile = sys.argv[1]
indexBaseDirectory = sys.argv[2]
outPrefix = sys.argv[3]
runDirectory = sys.argv[4]
queue = sys.argv[5]
nodes = sys.argv[6]
procs = sys.argv[7]
walltime = sys.argv[8]

# assemble cluster job
output, input = popen2('qsub')
input.write('#PBS -S /bin/bash\n')
input.write('#PBS -N ' + assembliesFile + '\n')
input.write('#PBS -d ' + runDirectory +'\n')
input.write('#PBS -e ' + runDirectory + ' -o ' + runDirectory + '\n')
input.write('#PBS -q ' + queue + '\n')
input.write('#PBS -l nodes=' + nodes + ':ppn=' + procs + ',walltime=' + walltime + '\n')

print("assembling transcripts from " + assembliesFile + " to " + outPrefix)

# commands for cluster job
input.write('cuffmerge ')
input.write('-g ' + indexBaseDirectory + '/Annotation/Genes/genes.gtf ')
input.write('-s ' + indexBaseDirectory + '/Sequence/Chromosomes ')
input.write('-p ' + procs + ' ') 
input.write('-o ' + outPrefix + ' ') 
input.write(assembliesFile + '\n')

#finish script
input.close()

