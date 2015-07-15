#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Convert any mixed fastq to 4 fastq files

import sys
import os

def read_read(infile, read_length):
    read = infile.readline()

    if not read or read == "":
        return None, None, None

    id1 = (read.split('/', 1)[0])[1:]
    strand = (read.split('/', 1)[1])[0]

    if strand != '2' and strand != '1':
        print("Wrong strands\n");
        return None, None, None

    delim = (read.split('/', 1)[0])[0]

    line = infile.readline()
    i = 0
    while line and i < (read_length - 1):
        read += line
        i += 1
        line = infile.readline()

    if not line:
        return id1, strand, read

    infile.seek(infile.tell() - len(line))
	
    return id1, strand, read
	


if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] + " <source>")	
    sys.exit()

inFileName = sys.argv[1]
inFile = open(inFileName, "r")

fName, ext = os.path.splitext(inFileName)
if ext == '.fq' or ext == '.fastq':
    read_length = 4
elif ext == '.fa' or ext == '.fasta':
    read_length = 2
else:
    print("Error: can't determine file format from extension! Please rename file to .fa, .fasta, .fq or .fastq")
    sys.exit()

pFile1 = open(fName + "_left" + ext, "w") 
pFile2 = open(fName + "_right" + ext, "w")
sFile1 = open(fName + "_single_left" + ext, "w") 
sFile2 = open(fName + "_single_right" + ext, "w")

id1, st1, read1 = read_read(inFile, read_length)
reads = {}

while id1 is not None:
	if id1 in reads:
		if st1 == '1':
			if reads[id1][1] != '2':
				print("Same strands in pair: " + str(id1) + "\n")
			
			pFile1.write(read1)
			pFile2.write(reads[id1][0])
		else:
			if reads[id1][1] != '1':
				print("Same strands in pair: " + str(id1) + "\n")
			
			pFile2.write(read1)
			pFile1.write(reads[id1][0])

		del reads[id1]
		
	else:
		reads[id1] = (read1, st1)

	id1, st1, read1 = read_read(inFile, read_length)


for id1, read in reads.iteritems():
	if read[1] == '1':
		sFile1.write(read[0])
	else:
		sFile2.write(read[0])


inFile.close()
pFile1.close()
pFile2.close()
sFile1.close()
sFile2.close()

