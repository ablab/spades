#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Generating insert size distribution from raw file

import sys

if len(sys.argv) < 5:
	print("Usage: <input raw file> <output> <min is> <max is> [fr/rf], fr -- default, rf -- calculate insert size for rf pairs as for chimeric")
	exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
minLen = int(sys.argv[3])
maxLen = int(sys.argv[4])

fr = True
if len(sys.argv) > 5 and sys.argv[5] == "rf":
	rf = False

hist = [0 for i in range(maxLen	+ 1)]

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
	else:
		cord = pos1 + len1 - pos2

	if cord < maxLen and cord > minLen:
		hist[cord] += 1

sum = 0
for i in range(minLen,maxLen):
	sum += hist[i]
        outFile.write(str(i) + ' ' + str(hist[i]) + '\n')

print("Total mate-pairs: " + str(sum))

inFile.close()
outFile.close()
