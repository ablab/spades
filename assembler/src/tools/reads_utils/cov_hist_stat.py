#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
maxLen = int(sys.argv[3])

hist = [0 for i in range(maxLen	+ 1)]

for line in inFile:
	stpos = int(line.split()[0])
	for i in range(0,int(line.split()[1])):
		cpos = stpos + i
		if cpos <= maxLen:
			hist[cpos] += 1

covered = 0.0
for i in range(0,maxLen + 1):
	if (hist[i] > 0):
		covered += 1.0
	outFile.write(str(i) + ' ' + str(hist[i]) + '\n')

print(str(covered/maxLen))

inFile.close()
outFile.close()
