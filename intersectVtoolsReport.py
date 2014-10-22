import argparse
import re
import pandas as pd

# add columns - unique-to-list1, intersect, unique-to-list2
# column values - sampleA:AF, sampleB:AF 
# Method:
# create dataframe for file A, and fileB
# melt dataframe 
# merge dataframe to create intersect, fileA only, and file B only
# Kyle Chang



def createDataFrame(filePath): 
	df = pd.read_csv(filePath, sep='\t')
	# remove columns after sample geno info
	df = df.iloc[0:,0:df.columns.get_loc('mutated_cases')]	
	# get a list of column names from chr to format
	ids=list(df.columns.values[:df.columns.get_loc('format')+1])
	vars=list(df.columns.values[df.columns.get_loc('format')+1:])
	# melt data frame - each row is a sample mutation
	df_melt = pd.melt(df, id_vars=ids, value_vars=vars, var_name='sample', value_name='info')
	# remove sample that doesn't have mutation - i.e. remove rows like FAP-EB1-P NA:::::
	df_melt = df_melt[~df_melt['info'].str.contains('NA')]
	return df_melt

def addAfColumn2DataFrame(df, info_name, format_name, var):
	# args: dataframe object, info column name, format column name, name of tool (mutect/torrent)

	adjaf = 'adj_allele_freq_' + var
	af = 'allele_freq_' + var
	# reset row index from 0 to n or adding new column will be affected by different row index
	df = df.reset_index(drop=True)
	# when merging dataframes, nan values need to be stored as float because pandas doesn't support int nan, thus chr and coord changed from int to float and it needs to be reset so it will print as int in csv. This is pandas gotchas.
	df['chr'] = df['chr'].astype(int)
	df['hg19_pos'] = df['hg19_pos'].astype(int)

	# splitting info string into a temp dataframe
	tmp = pd.DataFrame(df[info_name].str.split(':').tolist(), columns=list(df[format_name])[0].split(':'))
	# extract allele fraction from the last column and add to dataframe
	if 'AF' in tmp.columns:
		df[af] = tmp.loc[:,['AF']]
	elif 'FA' in tmp.columns:
		df[af] = tmp.loc[:,['FA']]
	else:
		raise Exception('missing allele fraction')
	
	df[adjaf] = tmp.loc[:,['AAF']]
	return df

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', dest='inputA')
	parser.add_argument('-b', dest='inputB')
	parser.add_argument('-aname', dest='nameA')
	parser.add_argument('-bname', dest='nameB')
	args = parser.parse_args()
	
	df_a = createDataFrame(args.inputA)
	df_b = createDataFrame(args.inputB)
	
	print 'file A only...'
	# union of file a and b
	union = pd.merge(df_a, df_b, on=['chr','hg19_pos','ref','alt','sample'], how='outer')
	# keep file a and columns only
	df_a_uniq = union[pd.notnull(union['info_x']) & pd.isnull(union['info_y'])].loc[:,'chr':'info_x']
	df_a_uniq = addAfColumn2DataFrame(df_a_uniq, 'info_x', 'format_x', args.nameA)

	print 'intersect...'
	# a and b overlap
	df_a_b = union[pd.notnull(union['info_x']) & pd.notnull(union['info_y'])]
	df_a_b = pd.concat([df_a_b.loc[:,'chr':'info_x'], df_a_b.loc[:,['format_y','info_y']]], axis=1)
	df_a_b = addAfColumn2DataFrame(df_a_b, 'info_x', 'format_x', args.nameA)
	df_a_b = addAfColumn2DataFrame(df_a_b, 'info_y', 'format_y', args.nameB)

	print 'file B only...'
	# union again, but list file b first, it's easier than picking out b columns from the first union
	union = pd.merge(df_b, df_a, on=['chr','hg19_pos','ref','alt','sample'], how='outer')
	df_b_uniq = union[pd.notnull(union['info_x']) & pd.isnull(union['info_y'])].loc[:,'chr':'info_x']
	df_b_uniq = addAfColumn2DataFrame(df_b_uniq, 'info_x', 'format_x', args.nameB)

	print 'write file...'
	# write files
	df_a_uniq.to_csv('fap_ampliseq.torrent-only.txt', sep='\t', index=False)
	df_a_b.to_csv('fap_ampliseq.intersect.txt', sep='\t', index=False)
	df_b_uniq.to_csv('fap_ampliseq.mutect-only.txt', sep='\t', index=False)

if __name__ == '__main__':
	main()
