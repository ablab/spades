#!/usr/bin/python -O

import sys

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
