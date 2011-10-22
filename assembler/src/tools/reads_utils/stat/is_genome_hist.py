#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
maxLen = int(sys.argv[3])
bar = int(sys.argv[4])

hist = {0:0}
rcount = {0:0}
for i in range(1,maxLen):
	hist[i] = 0
	rcount[i] = 0
	
for line in inFile:
	pos = int(line.split()[0]);
	hist[pos] += int(line.split()[1])
	rcount[pos] += 1

newHist = {0:0}
newRC = {0:0}
for i in range(0,maxLen/bar):
	newHist[i] = 0
	newRC[i] = 0

for i in range(0,maxLen):
	newHist[int(i/bar)] += hist[i]
	newRC[int(i/bar)] += rcount[i]

for i in range(0,maxLen/bar):
	if (newRC[i] != 0):
		outFile.write(str(i) + ' ' + str(newHist[i] / newRC[i]) + '\n')

inFile.close()
outFile.close()
