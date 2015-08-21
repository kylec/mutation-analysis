import re
import argparse

# function:
# - filter for de gene
# - split genes 
# removed # - flip log ratio sign (change from normal/polyp to polyp/normal)

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='infile')
parser.add_argument('-fs', dest='flip_sign', action='store_true')
parser.add_argument('-de', dest='de_switch', action='store_true')
parser.add_argument('-fc', dest='fc', type=float, default=1.0)
parser.add_argument('-qval', dest='qval', type=float, default=0.05)
args = parser.parse_args()

with open (args.infile, 'r') as genefile:
	# header
	print genefile.readline().rstrip()
	for line in genefile:
		cols = line.strip().split('\t')

		genes = cols[2].split(',')
		
		# if flip_sign switch is on
		if args.flip_sign:
			if re.match('-', cols[9]):
				cols[9] = cols[9].replace('-','')
			else:
				cols[9] = '-' + cols[9]

		# DE switch is on, filter for DE genes
		if args.de_switch:
			if float(cols[12]) > args.qval or (-args.fc < float(cols[9]) < args.fc): 
				continue
		
		# split genes
		for gene in genes:
			cols[2] = gene
			print '\t'.join(cols)

genefile.close()
		
