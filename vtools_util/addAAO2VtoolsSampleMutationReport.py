import argparse
import pandas as pd

# take a transformed vtools report (each row is a sample mutation) and add adjusted allele counts from adjusted allele frequency using DP
# Kyle Chang



def readDataFrame(filePath): 
	df = pd.read_csv(filePath, sep='\t')

	# splitting info string into a temp dataframe
	if 'info_y' in df.columns:
		info_name = 'info_y'
		format_name = 'format_y'
	else:
		info_name = 'info_x'
		format_name = 'format_x'

	tmp = pd.DataFrame(df[info_name].str.split(':').tolist(), columns=list(df[format_name])[0].split(':'))
	dfnew = pd.DataFrame()
	for row in tmp.iterrows():
		# iterrows returns a tuple , the data is actually in row[1]
		row = row[1]
		if ',' in row['AAF']:
			altCounts = str.split(row['AAF'],',')
			alt1 = round(float(altCounts[0])*float(row['DP']))
			alt2 = round(float(altCounts[1])*float(row['DP']))
			row['AAO'] = str(alt1) + ',' + str(alt2)
		else:
			alt = round(float(row['AAF'])*float(row['DP']))
			row['AAO'] = str(alt)
	
		dfnew = dfnew.append(row, ignore_index=True)

	df = pd.concat([df, dfnew], axis=1)
	return df

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', dest='inputA')
	args = parser.parse_args()
	
	df = readDataFrame(args.inputA)
	
	df.to_csv(args.inputA + '.aao.txt', sep='\t', index=False)

if __name__ == '__main__':
	main()
