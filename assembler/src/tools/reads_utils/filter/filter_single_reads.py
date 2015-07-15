#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Filter read files by bowtie log
#Leave only those that do not align

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
		return None, None

	infile.seek(infile.tell() - len(line))
	
	return id1, read
	


if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " <reads> <bowtie log>")	
	sys.exit()

rFileName1 = sys.argv[1]

rFile1 = open(rFileName1, "r")
blFile1 = open(sys.argv[2], "r")

ids = []

for line in blFile1:
	ids.append(line.split('/', 1)[0])

uids = set(ids)

fName1, ext1 = os.path.splitext(rFileName1)
outFile1 = open(fName1 + "_filtered" + ext1, "w") 

id1, read1 = read_read(rFile1)

while id1 is not None:
	if id1 not in uids:
		outFile1.write(read1)

	id1, read1 = read_read(rFile1)


rFile1.close()

outFile1.close()

blFile1.close()
