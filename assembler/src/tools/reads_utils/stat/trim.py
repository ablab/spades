#!/usr/bin/python -O

#Trim

import sys

if len(sys.argv) != 5:
	print("Usage: <input> <output> <lower limit> <upper limit>")
	exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
ll = int(sys.argv[3])
ul = int(sys.argv[4])

for line in inFile:
	coord = line.split(' ')
	x = int(coord[0])
	y = int(coord[1])
	
	if (x >= ll and x <= ul):
		outFile.write(str(x) + ' ' + str(y) + '\n' )

inFile.close()
outFile.close()
