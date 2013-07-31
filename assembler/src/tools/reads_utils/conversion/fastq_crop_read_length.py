#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Just for testing SPAdes on different read lengths: Crops FASTQ reads to MAX_LENGTH (all reads will be <= MAX_LENGTH)

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
if len(sys.argv) != 4:
	print("Usage: " + sys.argv[0] + " <input fastq> <output fastq> <max read length>")	
	sys.exit()

in_file = open(sys.argv[1])
out_file = open(sys.argv[2], 'w')
max_RL = int(sys.argv[3])

read = read_read(in_file)  # read[0] -- id, read[1] -- nucl. string, read[2] -- +, read[3] -- qual. string
while read:
    new_nucl_str = read[1]
    new_qual_str = read[3]
    if len(new_nucl_str) > max_RL:
        new_nucl_str = new_nucl_str[:max_RL] + '\n'
        new_qual_str = new_qual_str[:max_RL] + '\n'
    out_file.write(read[0])
    out_file.write(new_nucl_str)
    out_file.write(read[2])
    out_file.write(new_qual_str)
    read = read_read(in_file)

in_file.close()
out_file.close()
