#!/usr/bin/python
import sys
from  subprocess import call

inputFileName = sys.argv[1]
runDirectory = sys.argv[2]
procs = sys.argv[3]


print('Processing ' + inputFileName)

cmd = 'cufflinks -p ' + procs + ' -o ' + runDirectory + ' ' + inputFileName
print cmd
call(cmd, shell=True)
