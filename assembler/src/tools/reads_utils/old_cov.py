#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

delim = '/'

prevLine = ""
for line in inFile:
	if (prevLine == ""):
		prevLine = line
		continue

	if (line.split(delim, 1)[0] == prevLine.split(delim, 1)[0]):
		pos2 = int(line.split('\t', 4)[3])
		pos1 = int(prevLine.split('\t', 4)[3])
		outFile.write(str(pos1))
		outFile.write('\n')
                outFile.write(str(pos2))
                outFile.write('\n')
	else:
		print("Non-equal pairs\n")
		print(prevLine.split(delim, 1)[0])
		print(line.split(delim, 1)[0])

	prevLine = ""


inFile.close()
outFile.close()
