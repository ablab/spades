#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
chrsF = open(sys.argv[3])

chrs = {" ":0}
for line in chrsF:
        chrs[line.split()[0]] =	int(line.split()[1])

delim = '/'

for line in inFile:
	chr = line.split('\t',5)[2]
	pos1 = int(line.split('\t', 5)[3])
     	len1 = len(line.split('\t', 5)[4])
		
	outFile.write(str(pos1 + chrs[chr]) + ' ' + str(len1))
	outFile.write('\n')


inFile.close()
outFile.close()
