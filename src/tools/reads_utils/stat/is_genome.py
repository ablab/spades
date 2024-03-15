#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate insert size avergae value along the genome

import sys

if len(sys.argv) < 5:
	print("Usage: <input raw file> <output> <genome length> <histogram bar size> [fr/rf], fr -- default, rf -- calculate insert size for rf pairs as for chimeric")
	exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
maxLen = int(sys.argv[3])
bar = int(sys.argv[4])

fr = True
if len(sys.argv) > 5 and sys.argv[5] == "rf":
	rf = False

hist = [0 for i in range(maxLen	+ 1)]
rcount = [0 for i in range(maxLen + 1)]

while (1):
        line = inFile.readline()

        if not line:
                break

        pos1 = int(line.split(' ')[0])
        len1 = int(line.split(' ')[1])

        line = inFile.readline()

        if not line:
                break

        pos2 = int(line.split(' ')[0])
        len2 = int(line.split(' ')[1])

	if fr:
		cord = pos2 - pos1 + len2
		pos = (pos1 + pos2 + len2) / 2
	else:
		cord = pos1 + len1 - pos2
		pos = (pos1 + pos2 - len1) / 2

	hist[pos] += cord
	rcount[pos] += 1


newHist = [0 for i in range((maxLen + 1) / bar + 2)]
newRC = [0 for i in range((maxLen + 1) / bar + 2)]

for i in range(0,maxLen):
	newHist[int(i/bar)] += hist[i]
	newRC[int(i/bar)] += rcount[i]

for i in range(0,maxLen/bar):
	if (newRC[i] != 0):
		outFile.write(str(i) + ' ' + str(newHist[i] / newRC[i]) + '\n')

inFile.close()
outFile.close()
