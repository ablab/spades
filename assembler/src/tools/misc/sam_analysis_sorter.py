#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
import sys

if len(sys.argv) != 2:
    print ("SAM analyzer sorter")
    print ("Usage: " + sys.argv[0] + "  <input>")
    sys.exit(1)


inFile = open(sys.argv[1], "r")

hist = {}

maxLen = 0
line = inFile.readline()
while line:
	if not line.startswith("Contig") and not line.startswith("Total"):
		contigId = line.split()[0]
		hist[contigId] = int(line.split()[1])
		if len(contigId) > maxLen:
			maxLen = len(contigId)
	line = inFile.readline()

maxLen += 10

keys = hist.keys()
skeys = sorted(keys, key=lambda k: -hist[k])


title = "Contig ID"
spaces = ""
for i in range(0, maxLen - len(title)):
	spaces += " "
print (title + spaces + "Read number")

for contigId in skeys:
	spaces = ""
	for i in range(0, maxLen - len(contigId)):
		spaces += " "
	print (contigId + spaces + str(hist[contigId]))

inFile.close()


