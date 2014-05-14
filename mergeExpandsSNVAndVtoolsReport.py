import argparse

# Merge vtools annotation with expands subpopulation output
# Kyle Chang


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
#parser.add_arguemnt('-a', dest='anno')
parser.add_argument('-o', dest='output')
args = parser.parse_args()

# read annotation
# dictionary for annotation
anno=dict()
file = open('fap_mutect.report','r')
for line in file:
    data = line.rstrip('\n').split('\t')
    anno['\t'.join(data[0:2])] = '\t'.join([data[3],data[6], data[7],data[8],data[9]])

file.close()

# read subpopulation file
file = open(args.input,'r')
#output_file = open(args.output,'w')
#output_file.write(file.readline())
print file.readline().rstrip('\n') + "\t" + "\t".join(['dbSNP','region_type', 'region_name', 'mut_type', 'genename'])

for line in file:
    data = line.rstrip('\n').split('\t')
    key = '\t'.join(data[0:2])
    if key in anno:
        anno_info = anno[key]
    else:
        anno_info = ".\t.\t.\t.\t."
        
    #output_file.write('\t'.join([data,anno_info]))
    print line.rstrip('\n') + '\t' + anno_info 
    
file.close()
#output_file.close()