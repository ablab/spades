#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Trim Ns in fastq files by Hammer markers

import sys
import os
import shutil

def trim_read(infile):
	id1 = infile.readline()

	if not id1:
		return None

	seq = infile.readline()
	id2 = infile.readline()
	qual = infile.readline()

	if not seq or not id2 or not qual:
		print("Error in file")
		return None

	ltrim = 0
	lpos = id2.find("ltrim")
	if (lpos != -1):
		ltrim = int(id2[lpos + 6:].split()[0])

	rtrim = 0
	rpos = id2.find("rtrim")
	if (rpos != -1):
		rtrim = int(id2[rpos + 6:].split()[0])

	length = len(seq)

	if (len(seq) != len(qual)):
		print("Error in seq and quality lengths " + str(id1))

	return id1 + seq[ltrim:length - 1 - rtrim] + "\n" + id2 + qual[ltrim:length - 1 - rtrim] + "\n"

def trim_file(in_filename, out_filename):
    if in_filename == out_filename:    
        shutil.copy(in_filename, out_filename + "_temp")
        in_filename = out_filename + "_temp"

    inFile = open(in_filename, "r")
    outFile = open(out_filename, "w") 

    tread = trim_read(inFile)
    while tread is not None:
	    outFile.write(tread)
	    tread = trim_read(inFile)

    inFile.close()
    outFile.close()   
    
    if in_filename == out_filename + "_temp":
        os.remove(in_filename)

def main():
    if len(sys.argv) < 2:
	    print("Usage: " + sys.argv[0] + " <source>")	
	    sys.exit()

    inFileName = sys.argv[1]

    fName, ext = os.path.splitext(inFileName)
    outFileName = fName + "_trimmed" + ext 

    trim_file(inFileName, outFileName)

if __name__ == '__main__':
    main()
