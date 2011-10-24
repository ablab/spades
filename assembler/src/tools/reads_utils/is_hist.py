#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
maxLen = int(sys.argv[3])

hist = {0:0}
for i in range(1,maxLen):
	hist[i] = 0
	
for line in inFile:
	hist[int(line)] += 1

for i in range(0,maxLen):
	outFile.write(str(i) + ' ' + str(hist[i]) + '\n')

inFile.close()
outFile.close()
