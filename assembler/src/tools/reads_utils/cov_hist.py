#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
maxLen = int(sys.argv[3])
bar = int(sys.argv[4])


hist = [0 for i in range(maxLen + 1)]

for line in inFile:
	hist[int(line.split()[0])] = int(line.split()[1])	

newHist = [0 for i in range((maxLen + 1) / bar + 1)]

for i in range(maxLen + 1):
	newHist[int(i/bar)] += hist[i]

for i in range((maxLen + 1) / bar + 1):
	outFile.write(str(i) + ' ' + str(newHist[i] / bar) + '\n')

inFile.close()
outFile.close()
