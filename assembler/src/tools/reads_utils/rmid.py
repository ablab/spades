#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
file2 = open(sys.argv[2], 'w')


while 1:
	line = inFile.readline()

	if not line:
		break

	file2.write(line)

	line = inFile.readline()

	if not line:
		break

	file2.write(line)

	line = inFile.readline()

	if not line:
		break

	line = inFile.readline()

	if not line:
		break

	file2.write(line)
	
		
inFile.close()
file2.close()
