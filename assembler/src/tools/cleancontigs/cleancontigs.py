#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os

# Deletes reverse-complementary duplicates and simple duplicates from FASTA-file

sys.path.append(os.path.join(os.path.abspath(sys.path[0]), '..'))

import fastaparser

if len(sys.argv) < 2:
	print 'Usage', sys.argv[0], 'in.fasta > out.fasta'
	exit(1)

fastafilename = sys.argv[1]
fasta = fastaparser.read_fasta(fastafilename)
fasta_res = {}

for name, seq in fasta:
	if (seq not in fasta_res) and (fastaparser.rev_comp(seq) not in fasta_res):
		fasta_res[seq] = name

fastaparser.write_fasta((name, seq) for seq, name in fasta_res.iteritems())		
