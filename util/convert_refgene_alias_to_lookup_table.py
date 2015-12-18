import sys
from collections import defaultdict
from sets import Set

ref_gene_alias='/Users/kchang3/Analysis/references/refGene_alias.bed'
#ref_gene_alias='/Users/kchang3/Analysis/references/test/test_refGene_alias.bed'

# check if the gene has been read already
# check if alias is -

#print 'loading current symbols...'
current_genes = set()
with open (ref_gene_alias, 'r') as f:
    for line in f:
        cols = line.rstrip().split('\t')
        current_genes.add(cols[3])
        


#print 'mapping alias...'
# dictionary for mapping alias to current gene name
genes = dict()
with open (ref_gene_alias, 'r') as inputfile:
    for line in inputfile:
        cols = line.rstrip().split('\t')
        # 3 - gene, 4 alias
        current_name = cols[3]
            
        aliases = cols[4].split('|')

        for alias in aliases:
            # ignore if there is no alias 
            if alias == '-':
                next
            # ignore if alias is actually a current gene symbol, we trust the current symbol 
            elif alias in current_genes:
                #print 'already a gene symbol %s' % line
                next
            else:
                # map current gene to current gene
                genes[current_name] = current_name
                #check if alias been stored and see if they are pointing to the same gene
                if alias in genes.keys():
                    if genes[alias] != current_name:
                        #print 'WARNING: alias pointing to different gene name, %s' % alias 
                        genes[alias] = 'ambig'
                else:
                    genes[alias] = current_name
                    
#inputfile.closed
#print alias-gene dict
#print genes
print 'alias\tcurrent\n'
for key in genes.keys():
    print key + '\t' + genes[key] + '\n'