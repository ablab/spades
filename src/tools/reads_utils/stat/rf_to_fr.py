#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


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
