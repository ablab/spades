#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Filter read files by fastq
#Leave only those that do not appear in fastq file

import sys
import os

def read_readq(infile):
	read = infile.readline()

	if not read:
		return None, None, None

	id1 = read[1:]

	delim = read[0]

	#read
	read = infile.readline()
	#id
	infile.readline()
	#qual
	q = infile.readline()

	return id1, read, q


def read_readf(infile):
	read = infile.readline()

	if not read:
		return None, None

	id1 = read[1:]

	delim = read[0]

	#read
	read = infile.readline()

	return id1, read
	


if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " <sourse fastq> <corrected fastq>")	
	sys.exit()

rFileName1 = sys.argv[1]
rFileName2 = sys.argv[2]

rFile1 = open(rFileName1, "r")
rFile2 = open(rFileName2 , "r")

ids = {}

id1, read1, q = read_readq(rFile1)
while id1 is not None:
	ids[id1] = q
	id1, read1, q = read_readq(rFile1)


fName2, ext2 = os.path.splitext(rFileName2)
outFile = open(fName2 + ".fastq", "w") 

id1, read1 = read_readf(rFile2)
while id1 is not None:
	if id1 in ids:
		outFile.write("@" + id1 + read1 + "+" + id1 + ids[id1].strip()[0:len(read1) - 1])
		for i in range(len(ids[id1].strip()), len(read1) - 1):
			outFile.write("B") 
		outFile.write("\n")

	id1, read1 = read_readf(rFile2)


rFile1.close()
outFile.close()
rFile2.close()
