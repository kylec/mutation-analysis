import sys
from collections import defaultdict

#
ref_gene='/Users/kchang3/Analysis/references/refGene.bed'
gene_info='/Users/kchang3/Analysis/references/gene_info'
#ref_gene='/Users/kchang3/Analysis/references/testrefGene.bed'
#gene_info='/Users/kchang3/Analysis/references/testgene_info'

print "loading genes"
genes = defaultdict(set)
with open (ref_gene, 'r') as inputfile:
	for line in inputfile:
		cols=line.rstrip().split('\t')
		genes[cols[3]]

print "loading gene_info"
with open (gene_info, 'r') as inputfile:
	for line in inputfile:
		cols=line.rstrip().split('\t')
		if cols[2] in genes:
			for alias in cols[4].split('|'):
				genes[cols[2]].add(alias)

print "print genes"

with open (ref_gene, 'r') as inputfile, open(sys.argv[1], 'w') as outputfile:
    for line in inputfile:
		cols=line.rstrip().split('\t')
		alias = 'na'
		if cols[3] in genes:
			alias = '|'.join(genes[cols[3]]) 
		
		outputfile.write(line.rstrip() + '\t' + alias + '\n')
