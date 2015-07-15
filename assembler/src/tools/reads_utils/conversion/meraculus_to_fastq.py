#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Convert Meraculuos read file format to fastq

import sys

if len(sys.argv) != 3:
	print("Usage: <input> <output fastq> \n")
	exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

for line in inFile:
	if (line.strip() == ""):
		continue

	reads = line.strip().split(':')
	id1 = ""
	for i in range(0, len(reads) - 2):
		id1 += reads[i]

	outFile.write('@' + id1 + '\n')
	outFile.write(reads[len(reads) - 2] + '\n')
	outFile.write('+' + id1 + '\n')
	outFile.write(reads[len(reads) - 1] + '\n')

inFile.close()
outFile.close()
