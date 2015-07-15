#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import re

# Deletes reverse-complementary duplicates and simple duplicates from FASTA-file

sys.path.append(os.path.join(os.path.abspath(sys.path[0]), '..'))

import fastaparser

if len(sys.argv) < 2:
	print 'Usage:\n', sys.argv[0], '(-l|-c) (-n|-f) in.fasta > out.fasta'
	print "-l | -c\tsort by Length or Coverage"
	print "-n | -f\toutput Name only or Full recode"
	exit(1)

if sys.argv[1] == "-c":
	sortkey = lambda (name, seq): -float(re.search("cov(?:erage)?_([.0-9]+)", name).group(1))
else:
	sortkey = lambda (name, seq): -float(re.search("len(?:gth)?_([.0-9]+)", name).group(1))

nameonly = (sys.argv[2] == "-n")

fastafilename = sys.argv[3]
fasta = fastaparser.read_fasta(fastafilename)
fasta_res = []

for name, seq in fasta:
    fasta_res += [(name, seq)]

fasta_res = sorted(fasta_res, key=sortkey)

if nameonly:
    for name, seq in fasta_res:
    	print name
else:
	fastaparser.write_fasta((name, seq) for name, seq in fasta_res)
