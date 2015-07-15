#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Merge two single-read bowtie logs into paired one

import sys
import os

if len(sys.argv) != 3:
	print("Usage: " + sys.argv[0] + " <1st bowtie log> <2nd bowtie log>")	
	sys.exit()

rFileName1 = sys.argv[1]
rFileName2 = sys.argv[2]

rFile1 = open(rFileName1, "r")
rFile2 = open(rFileName2, "r")

pc = 0
c1 = 0
c2 = 0

ids = {}

for line in rFile1:
	id1 = line.split('/', 1)[0]
	ids[id1] = line
	c1 += 1
	
fName1, ext1 = os.path.splitext(rFileName1)
outFile = open(fName1 + "_paired" + ext1, "w")

for line in rFile2:
	c2 += 1
	id2 = line.split('/', 1)[0] 
	if id2 in ids:
		pc += 1
		outFile.write(ids[id2] + line)

print("1: " + str(c1) + " 2: " + str(c2) + " Pairs: " + str(pc) + "\n")


rFile1.close()
rFile2.close()

outFile.close()
