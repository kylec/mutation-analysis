#!/usr/bin/env python
#
# turn vtools mutect2 report into a sample mutation per line , add column "sample" and "af".
#
# usage:
# python2 mutect2ReportBySampleMut.py -i input.txt --exonic (keep exonic/splicing mutation only)

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
parser.add_argument('--exonic', dest='exonic',action='store_true')
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

  # skip non-exonic, splicing 
  if args.exonic:
    if fields[6] not in ['exonic','splicing']:
      continue
	
  mutfields = '\t'.join(fields[0:sampleIndex])

  # flag for reaching genotype column
  found = 0

  # ref,alt alleles
  refalt = fields[refIndex] + fields[altIndex]

  # index for matching sample name array and genotype info
  index = 0 
    
  # print samples
  for field in fields[sampleIndex:]:
    if field and not re.search('NA', field):
      dat = field.split(':')
      
      af = dat[3]
		
      print mutfields + '\t' + field + '\t' + samples[index] + '\t' + af
      
    index = index + 1

