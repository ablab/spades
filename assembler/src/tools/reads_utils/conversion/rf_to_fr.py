#!/usr/bin/python -O

import sys
import os
import itertools

def read_read(infile):
        read = infile.readline()

        if not read or read.strip() == "":
                return None, None, None

        id1 = (read.split('/', 1)[0])[1:]
        strand = (read.split('/', 1)[1])[0]

        if strand != '2' and strand != '1':
                print("Wrong strands\n");
                return None, None, None

        delim = (read.split('/', 1)[0])[0]

        line = infile.readline()
        i = 0
        while line and i < 3:
                read += line
                i += 1
                line = infile.readline()

        if not line:
                return id1, strand, read

        infile.seek(infile.tell() - len(line))
        
        return id1, strand, read


def comp(letter):
	return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}[letter.upper()]

def rev_comp(seq):
	return ''.join(itertools.imap(comp, seq[::-1]))

def rc_read(read):
	r = read.split('\n')
	return r[0] + "\n" + rev_comp(r[1])  + "\n" + r[2]  + "\n" + r[3][::-1]  + "\n"
        


if len(sys.argv) < 2:
        print("Usage: " + sys.argv[0] + " <source>")    
        sys.exit()


inFileName = sys.argv[1]
inFile = open(inFileName, "r")

fName, ext = os.path.splitext(inFileName)
outFile = open(fName + "_rc" + ext, "w") 

read = read_read(inFile)
while read[0]:
	outFile.write(rc_read(read[2]))
	read = read_read(inFile)

outFile.close()
inFile.close()
