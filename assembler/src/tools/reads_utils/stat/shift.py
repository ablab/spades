#!/usr/bin/python -O

#Shift both coordinates

import sys

if len(sys.argv) != 5:
	print("Usage: <input> <output> <x shift> <y shift>")
	exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')
shiftx = int(sys.argv[3])
shifty = int(sys.argv[4])

for line in inFile:
	coord = line.split(' ')
	x = int(coord[0])
	y = int(coord[1])
	
	outFile.write(str(shiftx + x) + ' ' + str(shifty + y) + '\n' )

inFile.close()
outFile.close()
