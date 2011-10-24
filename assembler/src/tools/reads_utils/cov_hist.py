#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
maxLen = int(sys.argv[3])
bar = int(sys.argv[4])

hist = {0:0}
for i in range(1,maxLen):
	hist[i] = 0
	
for line in inFile:
	for i in range(0,int(line.split()[1])):
		hist[int(line.split()[0]) + i] += 1

newHist = {0:0}
for i in range(0,maxLen/bar):
	newHist[i] = 0

for i in range(0,maxLen):
	newHist[int(i/bar)] += hist[i]

for i in range(0,maxLen/bar):
	outFile.write(str(i) + ' ' + str(newHist[i] / bar) + '\n')

inFile.close()
outFile.close()
