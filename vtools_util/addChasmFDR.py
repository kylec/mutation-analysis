import csv
import re
import sys

# Read chasm variant_result.tsv and add fdr value to vtools report 
# usage: python addChasmFDR.py [variant_result.tsv] [vtools_report.txt] [output file] [fdr|p]

colnum = 11 # chasm p-value column 0-based index
if sys.argv[3] == "fdr":
	colnum = 12

# read chasm
chasm = {}
with open(sys.argv[1], 'r') as chasmfile:
	for line in chasmfile:
		if re.search('chr\w+', line):
			cols = line.split('\t')
			chasm[re.sub('chr','',cols[2])+cols[3]+cols[5]+cols[6]] = cols[colnum]
			#print re.sub('chr','',cols[2])+cols[3]+cols[5]+cols[6] + ',' + cols[12]

# read report
with open (sys.argv[2], 'r') as reportfile, open (sys.argv[3], 'w') as outfile:
	reader = csv.DictReader(reportfile, delimiter='\t')	
	header = reader.fieldnames
	writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=header, lineterminator="\n")
	writer.writeheader()

	for row in reader:
		key = row['chr'] + row['hg19_pos'] + row['ref'] + row['alt']
		#print key
		if key in chasm:
			#print "found key"
			row['chasm_emprical_pval'] = chasm[key]

		writer.writerow(row)
	
