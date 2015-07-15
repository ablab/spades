#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import shutil

########################################################################	
### CONFIG & CHECKS 
########################################################################

fasta_width = 60
min_contig = 200
padding_length = 2000

import fastaparser

########################################################################

if len(sys.argv) < 3:
	print 'Contigs concatenator: makes one big contig from the assembly'
	print 'Usage: ', sys.argv[0], ' INPUT_FILE OUTPUT_FILE [COORDS_FILE]'
	sys.exit(0)

infilename = sys.argv[1]
outfilename = sys.argv[2]
coords = 0
if (len(sys.argv) > 3):
    coords = open(sys.argv[3], 'w')

padding = ""
for i in xrange(0, padding_length):
    padding += "N"

fasta = fastaparser.read_fasta(infilename)
summary_seq = ""
cur_coord = 1
for name, seq in fasta:
    if (len(seq) >= min_contig):
        if (coords != 0):
            coords.write(str(cur_coord) + " " + str(cur_coord + len(seq) - 1) + "\n")
            cur_coord += len(seq) + padding_length
        summary_seq += (seq + padding)

out = open(outfilename, 'w')
out.write(">sum_contig total_length=" + str(len(summary_seq)) + '\n')
for i in xrange(0,len(summary_seq),60):
    out.write(summary_seq[i:i+60] + '\n')

out.close()
if (coords != 0):
    coords.close()

print 'Done'
