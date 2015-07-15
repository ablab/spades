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

def read_read(infile):
	read = infile.readline()

	if not read:
		return None, None

	id1 = read[1:]

	delim = read[0]

	line = infile.readline()
	i = 0
       	while line and i < 3:
               	read += line
		i += 1
		line = infile.readline()

	if not line:
		return None, None

	infile.seek(infile.tell() - len(line))
	
	return id1, read
	


if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " <reads> <filtering fastq>")	
	sys.exit()

rFileName1 = sys.argv[1]

rFile1 = open(rFileName1, "r")
rFile2 = open(sys.argv[2], "r")

ids = []

id1, read1 = read_read(rFile2)
while id1 is not None:
	ids.append(id1)
	id1, read1 = read_read(rFile2)

uids = set(ids)

fName1, ext1 = os.path.splitext(rFileName1)
outFile1 = open(fName1 + "_filtered" + ext1, "w") 
outFile2 = open(fName1 + "_found" + ext1, "w")


id1, read1 = read_read(rFile1)
while id1 is not None:
	if id1 not in uids:
		outFile1.write(read1)
	else:
		outFile2.write(read1)

	id1, read1 = read_read(rFile1)


rFile1.close()
outFile1.close()
outFile2.close()
rFile2.close()
