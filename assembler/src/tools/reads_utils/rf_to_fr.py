#!/usr/bin/python -O

# Excahges read order by switching lines in raw file

import sys

if len(sys.argv) != 3:
	print("Usage: <raw file> <output>")
	exit(0)

inFile = open(sys.argv[1])
outFile = open(sys.argv[2], 'w')

while (1):
        line = inFile.readline()

        if not line:
                break

        line2 = inFile.readline()

        if not line2:
                break

	outFile.write(line2 + line)


inFile.close()
outFile.close()
