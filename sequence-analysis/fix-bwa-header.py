import re
import sys

# fix expanded tab in @PG group

with open (sys.argv[1], 'r') as f:
	for line in f:
		if re.match('@PG', line):
			lines = line.split('\t')
			#print lines
			#print lines[-1]
			line = '\t'.join(lines[0:5]) + '\\t' + '\\t'.join(lines[-5:])
			#print line	
		
		print line.rstrip('\n')
