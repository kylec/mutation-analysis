# sometimes a haploh output will have different chr start and end
# awk '{FS=OFS="\t"; print $9,$10,$13,$17,$15,$16}' eventcalls_postfilter_blood_as_ref.txt > eventcalls_postfilter_blood_as_ref.bed
# i.e. chr5 10000 chr6 100000
# result chr5 10000 end
#        chr6 start 10000


import csv
import sys
import copy 

# read chr length input
with open("/Users/kchang3/Analysis/references/chromInfo_nochr.txt", 'r') as inputfile:
	chrLen = dict(line.rstrip().split('\t') for line in inputfile)

#print chrLen['1']
#exit() 

with open(sys.argv[1], 'r') as csvfile, open(sys.argv[2], 'w') as outfile:
	reader = csv.DictReader(csvfile, delimiter='\t')
	header = reader.fieldnames
	writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=header)
	writer.writeheader()

	for cols in reader:
		if cols['start_chr'] != cols['end_chr']:
			#chrEnd = chrLen[cols['start_chr']]
			dupcols = copy.deepcopy(cols)
			dupcols['start_chr'] = cols['end_chr']
			dupcols['end_chr'] = cols['end_chr']
			dupcols['start_bp'] = 1
			dupcols['end_bp'] = cols['end_bp']
			# modify original chr
			cols['end_chr'] = cols['start_chr']
			# end bp is the last chr base of chrLen
			cols['end_bp'] = chrLen[cols['start_chr']]
			#print cols
			writer.writerow(cols)
			writer.writerow(dupcols)	
		else:
			writer.writerow(cols)
			#print cols

