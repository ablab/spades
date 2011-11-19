#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

delim = '/'

for line in inFile:
	pos1 = int(line.split('\t', 5)[3])
     	len1 = len(line.split('\t', 5)[4])
		
	outFile.write(str(pos1) + ' ' + str(len1))
	outFile.write('\n')


inFile.close()
outFile.close()
