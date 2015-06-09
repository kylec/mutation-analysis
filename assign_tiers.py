import csv
import sys
import re

# usage: python assign_tiers [input.txt] [crc_gene.txt] [output.txt]

#(1) Driver mutation: frameshift insertions and deletions, stop gains, stop loss mutations and mutations located in splice sites that were classified as a driver mutation by CHASM (based on an empirical p-value <= 0.05 from CHASM) and seen in multiple (2 or more) other tumors in the COSMIC database; 
#(2) Damaging recurrent mutations: predicted to be damaging by 2 or more algorithms and  seen in multiple (2 or more) other tumors in the COSMIC database( cosmicmutant count > 1) or in genes list  
#(3) Damaging mutations: predicted to be damaging by 2 or more algorithms; 
#(4) Potentially damaging mutations: predicted to be damaging by 1 algorithm; 
#(5) Passenger mutations: the remaining point mutations that are not stop gains or stop losses.

def is_damaging(s, damagingStrings):
	if s is None or s == "NA":
		return 0
	for damagingString in damagingStrings:
		if damagingString in s.split(';'):
			return 1
	else:
		return 0

# load crc genes list
with open (sys.argv[2], 'r') as genefile:
	gene_list = [line.strip() for line in genefile]

genefile.close()

with open (sys.argv[1], 'r') as csvfile:
	reader = csv.DictReader(csvfile, delimiter='\t')
	header = reader.fieldnames
	header = ["tier"] + header
	outfile = open (sys.argv[3], 'w')
	writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=header)
	writer.writeheader()

	for row in reader:
		tiers = 0
		
		# strip white space for mutations without pvalue
		chasm_pval = float(row['chasm_emprical_pval']) if row['chasm_emprical_pval'].strip() else 1.0
		cosmic_count = int(row['CosmicCodingMuts.CNT']) if row['CosmicCodingMuts.CNT'] != '.' else 0
		
		# prediction class
		lrtPred = is_damaging(row['LRT_pred'], ['D'])
		#condelPred = is_damaging(row['condel_pred'], ['D'])
		polyphenPred = is_damaging(row['Polyphen2_HDIV_pred'], ['D','P'])
		mutationTasterPred = is_damaging(row['MutationTaster_pred'], ['A','D'])
		siftScore = float(row['SIFT_score']) if row['SIFT_score'] != '.' else 1.0
		siftPred = 1 if siftScore < .05 else 0

		# count D in prediction programs
		#damage_count = lrtPred + condelPred + polyphenPred + mutationTasterPred
		damage_count = lrtPred + siftPred + polyphenPred + mutationTasterPred
		#damage_count = lrtPred + siftPred  + mutationTasterPred
		#damage_count = lrtPred + mutationTasterPred (actual run, not CNOT3 paper description)

		# Tier 1: indel OR stopgain/loss OR splicing OR missense(chasm < 0.05)
		# Tier 2: >=2 damage_count AND cosmic_count > 1 AND CRC gene list
		# Tier 3: >=2 damage_count
		# Tier 4: 1 damage_count 
		# Tier 5: else
		if re.search('stop|frameshift', row['mut_type']) or row['region_type'] == 'splicing' or (row['mut_type'] == 'nonsynonymous SNV' and chasm_pval < .05): 
			tiers = 1
		elif damage_count > 1 and (row['genename'] in gene_list or cosmic_count > 1):
			tiers = 2
		elif damage_count > 1:
			tiers = 3
		elif damage_count == 1:
			tiers = 4
		else:
			tiers = 5
		
		row["tier"] = tiers
		writer.writerow(row, )

csvfile.close()
outfile.close() 
