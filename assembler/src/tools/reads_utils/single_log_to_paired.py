#!/usr/bin/python

import sys
import os

if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " <1st bowtie log> <2nd bowtie log>")	
	sys.exit()

rFileName1 = sys.argv[1]
rFileName2 = sys.argv[2]

rFile1 = open(rFileName1, "r")
rFile2 = open(rFileName2, "r")

ids = {}

for line in rFile1:
	id1 = line.split('/', 1)[0]
	ids[id1] = line
	
fName1, ext1 = os.path.splitext(rFileName1)
outFile = open(fName1 + "_paired" + ext1, "w")

for line in rFile2:
	id2 = line.split('/', 1)[0] 
	if id2 in ids:
		outFile.write(ids[id2] + line)


rFile1.close()
rFile2.close()

outFile.close()
