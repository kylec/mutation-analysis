#/usr/bin/python
import sys
from subprocess import call

assembliesFile = sys.argv[1]
indexBaseDirectory = sys.argv[2]
outPrefix = sys.argv[3]
procs = sys.argv[4]


print("assembling transcripts from " + assembliesFile + " to " + outPrefix)

cmd = 'cuffmerge -g ' + indexBaseDirectory + '/Annotation/Genes/genes.gtf -s ' + indexBaseDirectory + '/Sequence/Chromosomes -p ' + procs + ' -o ' + outPrefix + ' ' + assembliesFile

print cmd
call(cmd, shell=True)

