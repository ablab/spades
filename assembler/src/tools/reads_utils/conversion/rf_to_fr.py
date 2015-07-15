#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import itertools

def read_read(infile):
        read = infile.readline()

        if not read or read.strip() == "":
                return None, None

        id1 = (read.split('/', 1)[0])[1:]

        line = infile.readline()
        i = 0
        while line and i < 3:
                read += line
                i += 1
                line = infile.readline()

        if not line:
                return id1, read

        infile.seek(infile.tell() - len(line))
        
        return id1, read


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
	outFile.write(rc_read(read[1]))
	read = read_read(inFile)

outFile.close()
inFile.close()
