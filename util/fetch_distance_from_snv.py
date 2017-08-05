#!/usr/bin/env python
#
# print snv position and  exon boundaries where the snv is sitting in
#
# usage:
# python2 fetch_distance_from_snv.py -i input.txt -bp 120 -tx ENST00000234420

import argparse
import re
import glob
from subprocess import check_output

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
parser.add_argument('-bp', dest='basepair')
parser.add_argument('-tx', dest='tx_id')
args = parser.parse_args()

ref = '/rsrch2/ccp_rsch/kchang3/references/ucsc.hg19.fasta'

# load exon boundaries
exon_dict = dict()
exons = glob.glob('/rsrch2/ccp_rsch/kchang3/kchang3/for_ester/Exon*.csv')
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
    #ref_allele = data[2]
    #var_allele = data[3]
    #tx_id = data[4]
    
    # lookup exon boundaries covering the snv
    exon_start = exon_end = 0
    for start, end in exon_dict[args.tx_id]:
        #print "looping list...", start,end
        if start <= snv_start & snv_start <= end:
            exon_start = start
            exon_end = end
            break
    
    #print snv_start, exon_start, exon_end
    if exon_start == 0 or exon_end == 0:
        raise Exception('ERROR: no exon found -' + line)   

    # print exon start, exon end, snv start, difference from start, difference from end
    print "%s\t%d\t%d\t%d\t%d\t%d" % (chrom, snv_start, exon_start, exon_end, snv_start-exon_start, exon_end - snv_start)
