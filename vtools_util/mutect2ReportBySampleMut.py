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
refIndex = samples.index('ref')
altIndex = samples.index('alt')
formatIndex = samples.index('format')
sampleIndex = formatIndex + 1
samples = samples[formatIndex+1:]

# print header
print '\t'.join(header.split('\t')[0:sampleIndex]) + '\t' + 'format_field' +  '\t' + "sample" + "\t" + "af"

for line in file:
  fields = line.rstrip('\n').split('\t')

  # skip non-exonic and synonymous
  #if (fields[6] == 'exonic' and fields[8] != 'synonymous SNV') or fields[6] == 'splicing':
  if (fields[6] == 'exonic' and fields[8] != 'unknown') or fields[6] == 'splicing':
  #if (1==1):
    # fields up to format field
    mutfields = '\t'.join(fields[0:sampleIndex])

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
    for field in fields[sampleIndex:]:
      if field and not re.search('NA', field):
        dat = field.split(':')
       
	
        #GT:GQ:AD:AF:ALT_F1R2:REF_F1R2:FOXOG:QSS
        af = dat[3]
		
        # indelocator - GT:AD_geno:DP_geno:N_DP_geno:T_DP_geno:N_AC_geno:T_AC_geno:N_SC_geno:T_SC_geno		
        #if re.search('-', refalt):
        #  counts = dat[8].split(',')
        #  af = (float(counts[0]) + float(counts[1])) / (float(counts[0]) + float(counts[1]) + float(counts[2]) + float(counts[3])) 
        #  af = str(round(af,3))

        # simplify samples
        #print "index= %s" % (index) 
        print mutfields + '\t' + field + '\t' + samples[index] + '\t' + af
      
      index = index + 1

