import re
import sys

# get a list of used genes/unsed genes for up/down regulated
# calculate percentage of used genes/unsed genes

up_genes_unused_count = 0
down_genes_unused_count = 0
up_genes_unused = ''
down_genes_unused = ''
up_genes_count = 0
down_genes_count = 0
up_genes = ''
down_genes = ''

with open (sys.argv[1], 'r') as ciddfile:
	for line in ciddfile:
		m1 = re.search('up-regulated signature genes not measured in the cmap drug expression data \((\w+)\): (.+)', line)
		if m1 != None:
			up_genes_unused_count =  float(m1.group(1))
			up_genes_unused =  m1.group(2)

		
		m2 = re.search('down-regulated signature genes not measured in the cmap drug expression data \((\w+)\): (.+)', line)
		if m2 != None:
			down_genes_unused_count = float(m2.group(1))
			down_genes_unused = m2.group(2)
		
		m3 = re.search('up-regulated signature genes \((\w+)\): (.+)', line)
		if m3 != None:
			up_genes_count = float(m3.group(1))
			up_genes = m3.group(2)

		m4 = re.search('Down-regulated signature genes \((\w+)\): (.+)', line)
		if m4 != None:
			down_genes_count = float(m4.group(1))
			down_genes = m4.group(2)

#print up_genes_count
#print up_genes_unused_count
#print down_genes_count
#print down_genes_unused_count

print "up-regulated genes used precentage: %.3f" % ((up_genes_count-up_genes_unused_count)/up_genes_count)
print "down-regulated genes used precentage: %.3f" % ((down_genes_count-down_genes_unused_count)/down_genes_count)
