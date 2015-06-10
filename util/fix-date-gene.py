import re
import sys

# if a file with line starts with 3-Mar
# replace it with Mar3


with open (sys.argv[1], 'r') as inputfile:
	header = inputfile.readline()
	print header.rstrip()
	for line in inputfile:
		m = re.match('(\d+)-(\w+)\t', line)
		if m != None:
			cols = line.rstrip().split('\t')
			cols[0] = m.group(2) + m.group(1)
			print '\t'.join(cols)
		else:
			print line.rstrip()

	
