# fetch flanking sequence a list of snv sites
# print fasta format
# > ref_seq
# ATATTTTTCGAGT
# > alt_seq
# ATTTAACGAGAGT
#
# Kyle Chang
import argparse
from subprocess import check_output

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
parser.add_argument('-bp', dest='basepair')
args = parser.parse_args()

ref = '/usr/local/epi/home/kchang3/references/ucsc.hg19.fasta'

file = open(args.input, 'r')
for line in file:
    data = line.rstrip('\n').split('\t')
    flank1_start = int(data[1]) - int(args.basepair)
    chrom = data[0].replace('chr','')
    flank1_end = int(data[1]) - 1
    flank2_start = int(data[1]) + 1
    flank2_end = int(data[1]) + int(args.basepair)
    flank1_seq = check_output('samtools faidx ' + ref + ' ' + 'chr' + chrom + ':' + str(flank1_start) + '-' + str(flank1_end) + ' | sed \'1d\'', shell=True)
    flank2_seq = check_output('samtools faidx ' + ref + ' ' + 'chr' + chrom + ':' + str(flank2_start) + '-' + str(flank2_end) + ' | sed \'1d\'', shell=True)
    # add 'wildtype to beginning of list data
    print '>' + '_'.join(['wildtype'] + data)
    print flank1_seq + data[2] + flank2_seq
    print '>' + '_'.join(['mutant'] + data)
    print flank1_seq + data[3] + flank2_seq
