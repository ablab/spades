#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Print reads overlaps from bowtie log

import sys

if len(sys.argv) != 3:
        print("Usage: <bowtie log file> <output>");
        exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

delim = '/'

prevLine = ""
for line in inFile:
	if (prevLine == ""):
		prevLine = line
		continue

	if (line.split(delim, 1)[0] == prevLine.split(delim, 1)[0]):
		l1 = line.split('\t', 6)
		l2 = prevLine.split('\t', 6)

		chr2 = l1[2]
		chr1 = l2[2]
		r1 = l2[3]
		r2 = l1[3]
		pos2 = int(r2)
		pos1 = int(r1)
               	len2 = len(l1[4])
               	len1 = len(l2[4])
		q2 = l1[5]
		q1 = l2[5]

		ovl = pos1 + len1 - pos2

		if ovl > 0 and chr1 == chr2:
			outFile.write('Overlap by ' + str(ovl) + ':\n')
			outFile.write(r1[:ovl] + '\n')
	                outFile.write(r2[-ovl:] + '\n')
                        outFile.write(q1[:ovl] + '\n')
                        outFile.write(q2[-olv:] + '\n')

	else:
		print("Non-equal pairs\n")
		print(prevLine.split(delim, 1)[0])
		print(line.split(delim, 1)[0])
		exit(0)

	prevLine = ""


inFile.close()
outFile.close()
