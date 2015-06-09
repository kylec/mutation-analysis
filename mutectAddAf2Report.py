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

# open file and get column index of fields
file = open(args.input, 'r')
header = file.readline()
samples = header.rstrip('\n').split('\t')
index = samples.index('format')
refIndex = samples.index('ref')
altIndex = samples.index('alt')
samples = samples[index+1:]

# allele fraction cutoff
afcutoff = .02

print header.rstrip('\n') + '\t' + "mutated_cases" + "\t" + "#_samples_>_af_cutoff"

for line in file:
  fields = line.rstrip('\n').split('\t')

  # skip non-exonic and synonymous
  #if (fields[6] == 'exonic' and fields[8] != 'synonymous SNV') or fields[6] == 'splicing':
  if (1==1):
    # empty string for mutated cases
    mutCase = ''
    # flag for reaching genotype column
    found = 0

    # ref,alt alleles
    refalt = fields[refIndex] + fields[altIndex]

    # index for matching sample name array and genotype info
    index = 0 
    # count for sample with allele fraction > .05
    sample_count = 0
    
    #print samples
    for field in fields:
      # found genotype column
      if re.search("GT:GQ", field):
        found = 1
        continue
   
      if found == 1:
        if field and not re.search('NA', field):
          dat = field.split(':')
          af = ''
          #GT:GQ:AD:DP:FA:BQ:SS
          af = dat[4]

          # check if af is greater than threshold
          if round(float(af),2) >= afcutoff:
            sample_count += 1

          # simplify samples
          #print "index= %s" % (index) 
          mutCase = mutCase + samples[index] + ':' + af + ','
        index = index + 1

    # print line only if at least one sample's af is greater than threshold
    print line.rstrip('\n') + "\t" + mutCase + "\t" + str(sample_count) 
