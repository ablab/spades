#!/usr/bin/python

import sys
import os
import re

# Deletes reverse-complementary duplicates and simple duplicates from FASTA-file

sys.path.append('../quality/libs')

import fastaparser

if len(sys.argv) < 2:
	print 'Usage', sys.argv[0], 'in.fasta > out.fasta'
	exit(1)

fastafilename = sys.argv[1]
fasta = fastaparser.read_fasta(fastafilename)
fasta_res = []

for name, seq in fasta:
    fasta_res += [(name, seq)]

fasta_res = sorted(fasta_res, key=lambda (name, seq): -float(re.search("cov_([.0-9]+)", name).group(1)))

fastaparser.write_fasta((name, seq) for name, seq in fasta_res)
