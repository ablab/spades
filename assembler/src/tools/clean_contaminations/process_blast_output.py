#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Read the result of BLAST and collect statistics on most common organisms.

import sys
import os
from collections import Counter

def readFile(fName, hitNumber):
	infile = open(fName, "r")
	line = infile.readline()
	datasets = ({}, [])
	while line:
		if line.startswith("Sequences producing significant alignments"):
			infile.readline()
			readSequences(infile, datasets, hitNumber)
		line = infile.readline()
	infile.close()
	return datasets
		
def readSequences(infile, datasets, hitNumber):
	line = infile.readline()
	i = 0
	while line != "\n" and i < hitNumber:
		i += 1
		data = line.strip().split('\t')
		seq = data[0]
		seqData = seq.strip().split("|")
		seqId = seqData[1]
		seqName = seqData[2].strip().split('  ')[0]
		(seqDict, seqList) = datasets
                if seqId not in seqDict:
                      seqDict[seqId] = seqName.strip()
                seqList.append(seqId)
		line = infile.readline()

def printStatistics(dataset, resultsNum):
        (dataDict, dataList) = dataset
        c = Counter(dataList)
        results = c.most_common(resultsNum)
        for (seqId, count) in results:
                print seqId, '\t', dataDict[seqId], '\t', count

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: " + sys.argv[0] + " <blast output> [<number of hits to consider>] [<number of results to output>]")	
		sys.exit()
	blastFileName = sys.argv[1]
	hitsNum = 20
	resultsNum = 30
	if len(sys.argv) > 2:
		hitsNum = int(sys.argv[2])
		if len(sys.argv) > 3:
			resultsNum = int(sys.argv[3])
	datasets = readFile(blastFileName, hitsNum)
	printStatistics(datasets, resultsNum)


