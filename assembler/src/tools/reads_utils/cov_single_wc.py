#!/usr/bin/python -O

#Create raw file for single reads using chromosome map file
#Chromosome map file should contain choromosome name (as used in bowtie) and its start position in genome, one chromosome per line

import sys

if len(sys.argv) != 3 && len(sys.argv) != 2:
        print("Usage: <bowtie log file> <output> [chromosomes file]");
        exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

wChrs = (len(sys.argv) == 3)

chrs = {" ":0}
if wChrs:	 
	chrsF = open(sys.argv[3])
	for line in chrsF:
		chrs[line.split()[0]] =	int(line.split()[1])

delim = '/'

for line in inFile:
	chr1 = line.split('\t',5)[2]
	pos1 = int(line.split('\t', 5)[3])
     	len1 = len(line.split('\t', 5)[4])

	addPos = 0
	if wChrs:
		addPos = chrs[chr1] 	
		
	outFile.write(str(pos1 + addPos) + ' ' + str(len1) + '\n')

inFile.close()
outFile.close()
