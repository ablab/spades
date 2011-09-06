#!/usr/bin/python

import sys
import os

# Deletes reverse-complementary duplicates and simple duplicates from FASTA-file

sys.path.append('../quality/libs')

import fastaparser

if len(sys.argv) < 2:
	print 'Usage', sys.argv[0], 'FASTA'
	exit(1)

fastafilename = sys.argv[1]
fasta = fastaparser.read_fasta(fastafilename)
fasta_res = {}

for name, seq in fasta:
	if (seq not in fasta_res) and (fastaparser.rev_comp(seq) not in fasta_res):
		fasta_res[seq] = name

fastaparser.write_fasta((name, seq) for seq, name in fasta_res.iteritems())

		

