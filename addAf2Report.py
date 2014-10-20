import argparse
import re

# add column with mutated cases
# filter:
#   skip non-exonic and synonymous
#   at least one sample has af above threshold
# Kyle Chang


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
args = parser.parse_args()

file = open(args.input, 'r')
header = file.readline()
samples = header.rstrip('\n').split('\t')
index = samples.index('format')
samples = samples[index+1:]
#print samples

print header.rstrip('\n') + '\t' + "mutated_cases" + "\t" + "#_samples_>_af_cutoff"

for line in file:
  fields = line.rstrip('\n').split('\t')

  # skip non-exonic and synonymous
  if (fields[6] == 'exonic' and fields[8] != 'synonymous SNV') or fields[6] == 'splicing':
    # empty string for mutated cases
    mutCase = ''
    # flag for reaching genotype column
    found = 0
    # ref,alt alleles
    refalt = fields[4] + fields[5]
    # index for matching sample name array and genotype info
    index = 0 
    # count for sample with allele fraction > .05
    sample_count = 0
    for field in fields:
      # found genotype column
      if re.search("GT:GQ", field):
        found = 1
        continue
   
      if found == 1:
        if not re.search('NA', field):
          dat = field.split(':')
          af = ''
          # complex case - when there are 2 alt counts. 
          # GT:GQ:RO:AO:DP:AF:AAF 1:50.0:945:61:1524:0.0423208:0.036213
          #FIXME: indel always the 2nd count?
      
				  # if indel, split AAF and choose 2nd count
          if re.search('-', refalt) and re.search(',', dat[3]):
            af = dat[-1].split(',')[1]
          # it's snv part of snv+indel, choose 1st count to be alt
          elif re.search(',', dat[3]):
            af = dat[-1].split(',')[0]      
          else:
            # simple case - there's only 1 alt cont
            af = dat[-1]

          # check if af is greater than > .05
          if round(float(af),2) >= 0.05:
            sample_count += 1

          # simplify samples
          mutCase = mutCase + samples[index].split('-')[1] + ':' + af + ','
        
        index = index + 1

    # print line only if at least one sample's af is greater than threshold
    print line.rstrip('\n') + "\t" + mutCase + "\t" + str(sample_count) 
