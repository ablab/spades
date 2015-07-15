#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Reads file <f.fasta> with multiple fasta sequences and builds an index <f.index>: "id start_position end_position".

import sys, os

def getId(line):
        data = line.strip().split("|")
        seqId = data[3]
	return seqId

def writePositions(fout, id, start, end):
	fout.write(id + '\t' + str(start) + '\t' + str(end) + '\n')

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "Usage: %s <name of file to process>" % sys.argv[0]
		sys.exit()
        finName = sys.argv[1]
	fin = open(finName, "r")
        fName, ext = os.path.splitext(finName)
        foutName = fName + ".index"
	fout = open(foutName, "w")
        line = fin.readline()
	id = None
	start = None
	end = None
        while line:
		if line.startswith(">"):
			end = fin.tell() - len(line)
			if id is not None and start is not None:			
				writePositions(fout, id, start, end)
			id = getId(line)
			start = end
		line = fin.readline()
	end = fin.tell()
	if id and start:
		writePositions(fout, id, start, end)
	fout.close()
	fin.close()
