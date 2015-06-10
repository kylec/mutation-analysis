import sys
import re

# function:
# - filter for de gene
# - split genes 
# - flip log ratio sign (change from normal/polyp to polyp/normal)

with open (sys.argv[1], 'r') as genefile:
	# header
	print genefile.readline().rstrip()
	for line in genefile:
		cols = line.strip().split('\t')

		genes = cols[2].split(',')
		if re.match('-', cols[9]):
			cols[9] = cols[9].replace('-','')
		else:
			cols[9] = '-' + cols[9]

		# filter
		if cols[6] == 'OK' and float(cols[12]) < .05 and (float(cols[9]) >= 1 or float(cols[9]) <= -1): 
			for gene in genes:
				cols[2] = gene
				print '\t'.join(cols)

genefile.close()
		
