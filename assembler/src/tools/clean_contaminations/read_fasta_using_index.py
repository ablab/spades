#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Given sequence id, find that sequence in a big multifasta file using prebuilt index.

import sys, os

def getId(line):
        data = line.strip().split("|")
        seqId = data[3]
	return seqId

def writePositions(fout, id, start, end):
	fout.write(id + '\t' + str(start) + '\t' + str(end) + '\n')

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "Usage: %s <id of fasta> <name of file to read from>" % sys.argv[0]
		sys.exit()
	id = sys.argv[1]
        ffastaName = sys.argv[2]
        fName, ext = os.path.splitext(ffastaName)
        findexName = fName + ".index"
	findex = open(findexName, "r")
	start = None
	end = None
	for line in findex:
		if line.startswith(id):
			data = line.strip().split()
			if id == data[0]:
				start = int(data[1])
				end = int(data[2])
				break
	findex.close()
	if start is not None and end is not None:
	        fin = open(ffastaName, "r")
		fin.seek(start)
		fasta = fin.read(end-start)
		print fasta,
		fin.close()
