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
		pos2 = int(line.split('\t', 5)[3])
		pos1 = int(prevLine.split('\t', 5)[3])
		r1 = line.split('\t', 5)[4]
		r2 = prevLine.split('\t', 5)[4]
                q1 = line.split('\t', 5)[5]
                q2 = prevLine.split('\t', 5)[5]
               	len2 = len(r1)
               	len1 = len(r2)

		ovl = pos1 + len1 - pos2

		if ovl > 0:
		
			outFile.write('Overlap by ' + str(ovl) + ':\n')
			outFile.write(r1[-ovl:])
			outFile.write('\n')
	                outFile.write(r2[:ovl+1])
	       	        outFile.write('\n')
                        outFile.write(q1[-ovl:])
                        outFile.write('\n')
                        outFile.write(q2[:ovl+1])
                        outFile.write('\n')

	else:
		print("Non-equal pairs\n")
		print(prevLine.split(delim, 1)[0])
		print(line.split(delim, 1)[0])

	prevLine = ""


inFile.close()
outFile.close()
