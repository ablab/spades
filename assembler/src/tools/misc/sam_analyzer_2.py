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
    print ("SAM analyzer")
    print ("Usage: " + sys.argv[0] + "  <SAM file>")
    sys.exit(1)


inFile = open(sys.argv[1], "r")

hist = {}
maxLen = 0

line = inFile.readline()
while line.startswith("@"):
	if line.startswith("@SQ"):
		contigId = line.strip().split()[1].split(':')[1]
		hist[contigId] = 0
		if len(contigId) > maxLen:
			maxLen = len(contigId)

	line = inFile.readline()

maxLen += 10


t = 0
while line:
	contigId = line.split()[1]
	hist[contigId] += 1
	line = inFile.readline()
	t += 1


title = "Contig ID"
spaces = ""
for i in range(0, maxLen - len(title)):
	spaces += " "
print (title + spaces + "Read number")


keys = hist.keys()
keys.sort()

ttl = 0

for contigId in keys:
	spaces = ""
	for i in range(0, maxLen - len(contigId)):
		spaces += " "
	print (contigId + spaces + str(hist[contigId]))
	ttl += hist[contigId]

print ("Total reads: " + str(t) + ", checksum: "+ str(ttl))

inFile.close()


