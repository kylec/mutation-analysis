# fetch flanking sequence a list of snv sites
# print fasta format
# > ref_seq
# ATATTTTTCGAGT
# > alt_seq
# ATTTAACGAGAGT
#
# Kyle Chang
import argparse
import re
import glob
from subprocess import check_output

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
parser.add_argument('-bp', dest='basepair')
args = parser.parse_args()

ref = '/usr/local/epi/home/kchang3/references/ucsc.hg19.fasta'

# load exon boundaries
exon_dict = dict()
exons = glob.glob('Exon*.csv')
for exon in exons:
    exon_file = open(exon, 'r')
    m = re.search('ENST\d+', exon)
    key = m.group(0)
    exon_list = []
    for line in exon_file:
        start, end = line.rstrip('\n').split('\t')
        exon_list.append((int(start),int(end)))
    # save key and list
    exon_dict[key] = exon_list

file = open(args.input, 'r')
for line in file:
    data = line.rstrip('\n').split('\t')
    flank1_start = int(data[1]) - int(args.basepair)
    chrom = data[0].replace('chr','')
    snv_start = int(data[1])
    ref_allele = data[2]
    var_allele = data[3]
    tx_id = data[4]
    
    # lookup exon boundaries covering the snv
    exon_start = exon_end = 0
    for start, end in exon_dict[tx_id]:
        #print "looping list...", start,end
        if start <= snv_start & snv_start <= end:
            exon_start = start
            exon_end = end
            break
    
    #print snv_start, exon_start, exon_end
    if exon_start == 0 or exon_end == 0:
        raise Exception('ERROR: no exon found -' + line)   

    # get flanking bases of exon boundaries
    flank1_start = exon_start - int(args.basepair)
    flank1_end = snv_start - 1
    flank2_start = snv_start + 1
    flank2_end = exon_end + int(args.basepair)
    flank1_seq = check_output('samtools faidx ' + ref + ' ' + 'chr' + chrom + ':' + str(flank1_start) + '-' + str(flank1_end) + ' | sed \'1d\'', shell=True)
    flank2_seq = check_output('samtools faidx ' + ref + ' ' + 'chr' + chrom + ':' + str(flank2_start) + '-' + str(flank2_end) + ' | sed \'1d\'', shell=True)
    # add 'wildtype to beginning of list data
    print '>' + '_'.join(['wildtype'] + data)
    print flank1_seq + data[2] + flank2_seq
    print '>' + '_'.join(['mutant'] + data)
    print flank1_seq + data[3] + flank2_seq
