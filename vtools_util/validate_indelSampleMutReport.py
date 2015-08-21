import sys
import csv

label = 'validation_prog'

with open(sys.argv[1], 'r') as inputfile, open(sys.argv[2], 'w') as outputfile:
	reader = csv.DictReader(inputfile, delimiter='\t')
	header = reader.fieldnames
	header = header + [label]
	writer = csv.DictWriter(outputfile, delimiter='\t', fieldnames=header, lineterminator="\n")
	writer.writeheader()

	for row in reader:
		dat = row['format_field'].split(':')
		#print dat

		# normal counts
		counts = dat[7].split(',')
		#af = (float(counts[0]) + float(counts[1])) / (float(counts[0]) + float(counts[1]) + float(counts[2]) + float(counts[3]))
		#af = str(round(af,3))
		if (int(counts[0]) + int(counts[1])) > 0 and row['thousandGenomes.EUR_AF_INFO'] == '.':
			row[label] = "FAIL"
		else:
			row[label] = "PASS"
		
		writer.writerow(row)

	
