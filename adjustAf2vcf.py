import re
import argparse
import pysam

# def get bam reads
def getBamReads(chr, start, end, ref, alt, bam):
	# thresholds
	fraction_base = 0.1
	ref_count_adj = 0
	alt_count_adj = 0
	ref_count = 0
	alt_count = 0
	af = 0
	af_adj = 0

	# read bam
	samfile = pysam.Samfile(bam, 'rb')
	for pileupcolumn in samfile.pileup(chr, start, end):
		if pileupcolumn.pos == start:
			for pileupread in pileupcolumn.pileups:
				
				# extract pileup bases 
				# pileupread.qpos includes softclip, thus use pileupread.alignment.seq becauase it includes softclip 
				seq = pileupread.alignment.seq[pileupread.qpos]
				
				# pileupread.alignment.qlen is aligned length, does NOT include softclip
				#print '\tbase in read %s, pos=%s, base=%s, read_start=%s, len=%s' % (pileupread.alignment.qname, pileupread.qpos, seq, pileupread.alignment.pos, pileupread.alignment.qlen)
				exclude_5prime = float(pileupread.alignment.qlen) * fraction_base
				exclude_3prime = float(pileupread.alignment.qlen) - float(pileupread.alignment.qlen) * fraction_base 
				#print  '5prime=%s, 3prime=%s' % (exclude_5prime, exclude_3prime)
			
				# count ref and alt occurence
				if ref == seq: 
 					ref_count += 1
				elif alt == seq:
					alt_count += 1 
			
				# current pileup base index includes softclip	
				# should I include softclip when excluding variant base?
				if (pileupread.qpos > exclude_5prime) and (pileupread.qpos < exclude_3prime):
					if ref == seq:
						ref_count_adj += 1
					elif alt == seq:
						alt_count_adj += 1
				
	samfile.close()
  
	if ref_count + alt_count > 0:
		af = float(alt_count) / float(ref_count + alt_count)
	if ref_count + alt_count_adj > 0:	
		af_adj = float(alt_count_adj) / float(ref_count + alt_count_adj)
	
	#print 'site=%s:%s,%s/%s, before=%s,%s,%s\tafter=%s,%s,%s' % (chr, end, ref, alt, ref_count, alt_count, af, ref_count, alt_count_adj, af_adj)
	return af_adj
   
# main read vcf
# extract position
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', dest='vcf')
	parser.add_argument('-b', dest='bam')
	parser.add_argument('-s', dest='sample')
	args = parser.parse_args()
  
# do not use vcf reader. values can't be edited or follow here https://github.com/jamescasbon/PyVCF/issues/82	
	file = open(args.vcf, 'r')
	for line in file:
		if re.match('^#', line):
			print line.rstrip('\n')
		else:
			# convert to 0-base coordinates
			dat = line.rstrip('\n').split('\t')
			start = int(dat[1]) - 1
			end = int(dat[1])
		
			# skip multi_allelic mutations
			# does indel work?
			af = '.'
			if re.search(',', dat[4]) or len(dat[3]) > 1 or len(dat[4]) > 1:  
				af = '.'
			else:
				af = getBamReads(dat[0], start, end, dat[3], dat[4], args.bam)
		
			# format
			dat[8] = 'AAF:' + dat[8]
			dat[9] = str(af) + ':' + dat[9]
			print '\t'.join(dat)
	file.close()

if __name__ == '__main__':
	main()

