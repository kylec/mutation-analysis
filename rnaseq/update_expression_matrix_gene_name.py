import sys

# usage: update gene names in a expression matrix , assuming first column is the gene name
# python update_expression_matrix_gene_name.py [gene_alias_lookup.txt] [output.txt]

# read gene lookup table
# read expression matrix
# update gene if it's a key in alias table
# dict for alias and genes
genes = {}
gene_lookup_file = '/Users/kchang3/Analysis/references/gene_alias_lookup.txt'
with open (gene_lookup_file, 'r') as f:
    for line in f:
        (key, value) = line.rstrip().split('\t')
        genes[key] = value
        
#print genes

# read 
with open (sys.argv[1], 'r') as f, open(sys.argv[2], 'w') as outputfile:
    #print header
    outputfile.write(f.readline())
    for line in f:
        status = 'noupdate'
        cols = line.rstrip().split('\t')
        # if gene is an alias , get current name
        curr_name = 'na'
        if cols[0] in genes.keys():
            curr_name = genes[cols[0]]
            if genes[cols[0]] == "ambig":
                status = 'ambig'
            else:
                cols[0] = curr_name
                status = 'update'
                
        # log whether a line is updated
        print cols[0] + '\t' + curr_name + '\t' + status
            
        outputfile.write('\t'.join(cols) + '\n')
