import sys
from collections import defaultdict

# usage: add gene alias to refGene.bed
# example: python add_ref_genes_alias.py [output_file]

ref_gene='/Users/kchang3/Analysis/references/refGene.bed'
gene_info='/Users/kchang3/Analysis/references/gene_info'
#ref_gene='/Users/kchang3/Analysis/references/test/testrefGene.bed'
#gene_info='/Users/kchang3/Analysis/references/test/testgene_info'

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
#print genes

with open (ref_gene, 'r') as inputfile, open(sys.argv[1], 'w') as outputfile:
    for line in inputfile:
        cols=line.rstrip().split('\t')

        # if gene has alias
        if len(genes[cols[3]]) > 0:
            alias = '|'.join(genes[cols[3]]) 
        else:
            alias = '-'

        outputfile.write(line.rstrip() + '\t' + alias + '\n')
