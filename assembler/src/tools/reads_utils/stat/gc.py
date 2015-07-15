#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Generating GC distribution of genome and GC/insert size correlation from raw file

import sys
import os

def get_avg(arr, start, end):
	s = 0
	for i in range(start, end + 1):
		s += arr[i]

	return int (float(s) / float(end - start + 1))


if len(sys.argv) < 8:
	print("Usage: <input raw file> <genome> <genome length> <bar size> <min is> <max is> <--window|--gap|--reads|--pair> [fr/rf], fr -- default, rf -- calculate insert size for rf pairs as for chimeric")
	exit(0)

infName = sys.argv[1]
inFile = open(infName)
genome = open(sys.argv[2])
genomeLen = int(sys.argv[3])
bar = int(sys.argv[4])
minins = int(sys.argv[5])
maxins = int(sys.argv[6])

method = 0
if sys.argv[7] == "--reads":
	method = 1
elif sys.argv[7] == "--pair":
	method = 2
elif sys.argv[7] == "--gap":
	method = 3

if method != 0:
	bar = 1

fr = True
if len(sys.argv) > 8 and sys.argv[8] == "rf":
	rf = False

ACGT = "ACGT"
GC = "GC"

fName, ext = os.path.splitext(infName)
gcFile = open(fName + "_w" + str(bar) + ".gc", "w")
gcHist = [0.0 for i in range(genomeLen/bar + 1)]

genome.readline();
count = 0
gc = 0
pos = 0
line = genome.readline()
while line:
	for c in line:
		if c in ACGT:
			count += 1
			if c in GC:
				gc += 1

		if count == bar:
			gcHist[pos] = 100.0 * float(gc)/float(count)
			pos += 1
			gc = 0
			count = 0

	line = genome.readline()
if count != 0:
	gcHist[pos] = 100.0 * float(gc)/float(count)

for i in range(genomeLen/bar + 1):
	gcFile.write(str(i) + " " + str(gcHist[i]) + "\n")

gcFile.close()
genome.close()

maxpercent = 100
gcisFile = open(fName + "_w" + str(bar) + "_m" + str(method) + ".gcis", "w")
gcisHist = [[0 for i in range(maxins - minins + 1)] for j in range(maxpercent + 1)]

line = inFile.readline()
while line:
	l1 = line.split(' ')
	pos1 = int(l1[0])
	len1 = int(l1[1])

        line = inFile.readline()

        if not line:
                break

	l2 = line.split(' ')
       	pos2 = int(l2[0])
       	len2 = int(l2[1])


	if fr:
		inss = pos2 - pos1 + len2
		pos = pos1 + inss/2
	else:
		inss = pos1 + len1 - pos2
		pos = pos2 + inss/2

	gc_percent = 0
	if method == 0:
		gc_percent = int(gcHist[pos/bar])
	elif method == 1:
		gc_percent = int(100.0 * (float(l1[2]) + float(l2[2])) / float(len1 + len2))
	elif method == 2:
		gc_percent = get_avg(gcHist, pos1, pos2 + len2)
	elif method == 3:
		gc_percent = get_avg(gcHist, pos1 + len1, pos2)

	if inss < maxins and inss >= minins and pos < genomeLen:
		gcisHist[gc_percent][inss - minins] += 1

	line = inFile.readline()


for i in range(maxpercent + 1):
	for j in range(maxins - minins + 1):
	        gcisFile.write(str(j) + ' '  + str(i) + ' ' + str(gcisHist[i][j]) + '\n')
	

inFile.close()
gcisFile.close()
