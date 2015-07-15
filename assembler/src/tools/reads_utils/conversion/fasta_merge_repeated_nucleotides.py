#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Convert contigs (i.e a reference) for experiment of running SPAdes on E. coli MC reads in "IonTorrent" mode
# (all series of repeated nucleotides are changed to single nucleotides).

import sys
import os
import fastaparser

# MAIN
if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " <input fasta> <output fasta>")	
	sys.exit()

new_fasta = []
for name, seq in fastaparser.read_fasta(sys.argv[1]): 
    new_seq = seq[0]
    for i in range(1, len(seq)):
        if seq[i - 1] != seq[i]:
            new_seq += seq[i]
    new_fasta.append((name, new_seq))

fastaparser.write_fasta_to_file(sys.argv[2], new_fasta)
