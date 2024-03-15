#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Create raw file for paired reads using chromosome map file (optional)
#Chromosome map file should contain choromosome name (as used in bowtie) and its start position in genome, one chromosome per line

import sys


def count_gc(read):
	GC = "GC"
	gc = 0
	for c in read:
		if c in GC:
			gc += 1

	return gc

if len(sys.argv) != 4 and len(sys.argv) != 3:
        print("Usage: <bowtie log file> <output> [chromosomes file]");
        exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

wChrs = (len(sys.argv) == 4)

chrs = {" ":0}
if wChrs:	 
	chrsF = open(sys.argv[3])
	for line in chrsF:
		chrs[line.split()[0]] =	int(line.split()[1])
	chrsF.close()


delim = '/'

prevLine = ""
for line in inFile:
	if (prevLine == ""):
		prevLine = line
		continue

	if (line.split(delim, 1)[0] == prevLine.split(delim, 1)[0]):
		l1 = line.split('\t', 5)
		l2 = prevLine.split('\t', 5)

		chr2 = l1[2]
		chr1 = l2[2]
		pos2 = int(l1[3])
		pos1 = int(l2[3])
               	len2 = len(l1[4].strip())
               	len1 = len(l2[4].strip())

		if not wChrs or chr1 == chr2:
			addPos = 0
			if wChrs:
				addPos = chrs[chr1] 	

			outFile.write(str(pos1 + addPos) + ' ' + str(len1) + ' ' + str(count_gc(l2[4])))
			outFile.write('\n')
		        outFile.write(str(pos2 + addPos) + ' ' + str(len2) + ' ' + str(count_gc(l1[4])))
		        outFile.write('\n')
	else:
		print("Non-equal pairs\n")
		print(prevLine.split(delim, 1)[0])
		print(line.split(delim, 1)[0])
		exit(0)

	prevLine = ""

inFile.close()
outFile.close()
