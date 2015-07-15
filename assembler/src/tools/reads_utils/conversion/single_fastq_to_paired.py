#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Merge two fastq files into paired one

import sys
import os

def read_read(infile):
       	read = infile.readline()

       	if not read:
               	return None, None

       	id1 = (read.split('/', 1)[0])[1:]

       	delim = (read.split('/', 1)[0])[0]

       	line = infile.readline()
	i = 0
       	while line and i < 3:
               	read += line
		i += 1
               	line = infile.readline()

       	if not line:
               	return id1, read

       	infile.seek(infile.tell() - len(line))

       	return id1, read


if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " <1st fastq> <2nd fastq>")	
	sys.exit()

rFileName1 = sys.argv[1]
rFileName2 = sys.argv[2]

rFile1 = open(rFileName1, "r")
rFile2 = open(rFileName2, "r")

pc = 0
c1 = 0
c2 = 0

ids = {}

id1, read1 = read_read(rFile1)
while id1 is not None:
	ids[id1] = read1
	c1 += 1
	id1, read1 = read_read(rFile1)
	
fName1, ext1 = os.path.splitext(rFileName1)
outFile = open(fName1 + "_paired" + ext1, "w")

id2, read2 = read_read(rFile2)
while id2 is not None:
	c2 += 1
	if id2 in ids:
		pc += 1
		outFile.write(ids[id2] + read2)

	id2, read2 = read_read(rFile2)


print("1: " + str(c1) + " 2: " + str(c2) + " Pairs: " + str(pc) + "\n")


rFile1.close()
rFile2.close()

outFile.close()
