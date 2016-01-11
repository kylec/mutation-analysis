import sys
import csv

label = 'validation'

with open(sys.argv[1], 'r') as inputfile, open(sys.argv[2], 'w') as outputfile:
	reader = csv.DictReader(inputfile, delimiter='\t')
	header = reader.fieldnames
	if label not in header:
		header = header + [label]
	writer = csv.DictWriter(outputfile, delimiter='\t', fieldnames=header, lineterminator="\n")
	writer.writeheader()

	for row in reader:
		if row['validation_visual'] == "FAIL" or row['validation_prog'] == "FAIL": 	
			row[label] = "FAIL"
		else:
			row[label] = "PASS"	
			
		writer.writerow(row)

	
