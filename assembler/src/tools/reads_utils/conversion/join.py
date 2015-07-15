#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Merge two read files into one

import sys

if len(sys.argv) < 4:
	print("Usage: <1st file> <2nd file> <output> [format: fastq/fasta] \n")
	exit(0)

file1 = open(sys.argv[1])
file2 = open(sys.argv[2])
outFile = open(sys.argv[3], 'w')

format = 'fastq'
lines = 4
if len(sys.argv) == 5:
	format = sys.argv[4]

if format == 'fasta':
	lines = 2
	

while 1:
	line1 = file1.readline()
	line2 = file2.readline()

	if not line1:
		break
	if not line2:
		break

	if (line1.split('/', 1)[0] != line2.split('/', 1)[0]):
		print 'Error comaring indices'
		break

	outFile.write(line1)
	for i in range(1,lines):
		l = file1.readline()
			
		if not l:
			break
		outFile.write(l)

	outFile.write(line2)
	for i in range(1,lines):
		l = file2.readline()
			
		if not l:
			break
		outFile.write(l)
		
file1.close()
file2.close()
outFile.close()
