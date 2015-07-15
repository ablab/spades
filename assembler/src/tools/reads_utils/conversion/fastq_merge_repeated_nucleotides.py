#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Convert reads for experiment of running SPAdes on E. coli MC reads in "IonTorrent" mode
# (all series of repeated nucleotides changed to single nucleotides).

import sys
import os

def read_read(in_file):
    read = []
    read_id = in_file.readline()
    if not read_id:
       	return []

    read.append(read_id)
    for i in range(3):
        read.append(in_file.readline())    

    return read


# MAIN
if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " <input fastq> <output fastq>")	
	sys.exit()

in_file = open(sys.argv[1])
out_file = open(sys.argv[2], 'w')

read = read_read(in_file)  # read[0] -- id, read[1] -- nucl. string, read[2] -- +, read[3] -- qual. string
while read:
    new_nucl_str = read[1][0]
    new_qual_str = read[3][0]
    for i in range(1, len(read[1])):
        if read[1][i - 1] != read[1][i]:
            new_nucl_str += read[1][i]
            new_qual_str += read[3][i]
    out_file.write(read[0])
    out_file.write(new_nucl_str)
    out_file.write(read[2])
    out_file.write(new_qual_str)
    read = read_read(in_file)

in_file.close()
out_file.close()
