#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Convert any unordered 2 fastqs to 3 fastq files -- 2 paired and 1 single
# does not imply /1 or /2 in reads ids

import sys
import os

def read_read(infile, read_length):
    read = infile.readline()

    if not read or read == "":
        return None, None

    id1 = (read.split('/', 1)[0])[1:]

    line = infile.readline()
    i = 0
    while line and i < (read_length - 1):
        read += line
        i += 1
        line = infile.readline()

    if not line:
        return id1, read

    infile.seek(infile.tell() - len(line))
	
    return id1, read
	


if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <left reads> <right reads> <output>")	
    sys.exit()


leftFileName = sys.argv[1]
rightFileName = sys.argv[2]
leftFile = open(leftFileName, "r")
rightFile = open(rightFileName, "r")


fName, ext = os.path.splitext(leftFileName)
if ext == '.fq' or ext == '.fastq':
    read_length = 4
elif ext == '.fa' or ext == '.fasta':
    read_length = 2
else:
    print("Error: can't determine file format from extension! Please rename file to .fa, .fasta, .fq or .fastq")
    sys.exit()


outFileName = sys.argv[3]
pFile1 = open(outFileName + "_left" + ext, "w") 
pFile2 = open(outFileName + "_right" + ext, "w")
sFile = open(outFileName + "_single" + ext, "w") 


id1, read1 = read_read(leftFile, read_length)
left_reads = {}

while id1 is not None:
	left_reads[id1] = read1
	id1, read1 = read_read(leftFile, read_length)


id2, read2 = read_read(rightFile, read_length)
right_reads = []

while id2 is not None:
	if not id2 in left_reads:
		right_reads.append(read2)
	else:
		pFile1.write(left_reads[id2])
		pFile2.write(read2)
		del left_reads[id2]

	id2, read2 = read_read(rightFile, read_length)


for id1, read1 in left_reads.iteritems():
	sFile.write(read1)

for read2 in right_reads:
	sFile.write(read2)


leftFile.close()
rightFile.close()
pFile1.close()
pFile2.close()
sFile.close()
