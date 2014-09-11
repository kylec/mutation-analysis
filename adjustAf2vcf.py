import re
import argparse
import pysam

#NOTE
# pileupcolumn.n doesn't get you all the read counts because pileup has to iterate over # all columns, weird.

# def get bam reads
def getBamReads(chr, start, end, ref, alt, bam):
	# thresholds
	fraction_base = 0.1
	alt_count = 0
	read_count = 0
	af = 0.0

	# read bam
	samfile = pysam.Samfile(bam, 'rb')
	for pileupcolumn in samfile.pileup(chr, start, end):
		# reach mutation site
		if pileupcolumn.pos == start:
			for pileupread in pileupcolumn.pileups:
				# pileupread.qpos includes softclip, thus use pileupread.alignment.seq becauase it includes softclip 
				# length of allele
				baselen = max(len(ref), len(alt))
				# extract pileup bases 
				seq = pileupread.alignment.seq[pileupread.qpos:pileupread.qpos+baselen]
				#print 'args' + ref+alt + ','+ str(baselen) + ','+ pileupread.alignment.seq[pileupread.qpos] + ',' + seq + ',indel=' + str(indel_count) + ',' + pileupread.alignment.qname + ',' + pileupread.alignment.seq + ', ' + str(pileupread.qpos) + ',indtrue' + str(pileupread.indel) + ',' + str(start) + ',' + pileupread.alignment.cigarstring	
				
				# pileupread.alignment.qlen is aligned length, does NOT include softclip
				#NOTE should include softclip?
				exclude_5prime = float(pileupread.alignment.qlen) * fraction_base
				exclude_3prime = float(pileupread.alignment.qlen) - exclude_5prime
			
				# count alt if it's not at the edge
				if (pileupread.qpos > exclude_5prime) and (pileupread.qpos < exclude_3prime):
					if alt == seq and len(alt) == len(ref):
						# snv
 						alt_count += 1
					elif len(alt) > 1:
					# insertion - check allele match , pileupread.indel can be any positive int
						if alt == seq and pileupread.indel > 0:
							alt_count += 1
					elif len(ref) > 1:
						# deletion
						#FIXME ref!=seq may not be exactly the same deleted base, close enough? 
						if ref != seq and pileupread.indel < 0:
							alt_count += 1
				
				# count read
				read_count += 1
	
	samfile.close()
  
	af = round(float(alt_count) / float(read_count), 6)
	#print 'site=%s:%s,%s/%s, counts=%s,%s,%s' % (chr, end, ref, alt, read_count, alt_count, af)
	return str(af)
   
# main read vcf
# extract position
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', dest='vcf')
	parser.add_argument('-b', dest='bam')
	args = parser.parse_args()
  
# do not use vcf reader. values can't be edited or follow here https://github.com/jamescasbon/PyVCF/issues/82	
	file = open(args.vcf, 'r')
	for line in file:
		if re.match('^#', line):
			print line.rstrip('\n')
		else:
			# convert to 0-base coordinates
			dat = line.rstrip('\n').split('\t')
			chr = dat[0]
			start = int(dat[1]) - 1
			end = int(dat[1])
			ref = dat[3]
			alt = dat[4]

			# skip multi_allelic mutations
			# does indel work?
			af = '.'
			# skip multi-alt
			if alt.count(',') > 1:
				continue				
			elif alt.count(',') == 1:
				# 2 allele is ok
				# list of alt allele and af
				alts = alt.split(',')
				afs = []
				for alt in alts:
					#print 'args' + ref + alt
					af = getBamReads(chr, start, end, ref, alt, args.bam)
					# push af to a list, join by ,
					afs.append(str(af))
				af = ','.join(afs)
			else:
				# snv
				af = getBamReads(chr, start, end, ref, alt, args.bam)
		
			# format
			dat[8] = dat[8] + ':AAF' 
			dat[9] = dat[9] + ':' + af
			print '\t'.join(dat)
	file.close()

if __name__ == '__main__':
	main()

