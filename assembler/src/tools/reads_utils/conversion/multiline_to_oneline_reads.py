#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import itertools
import shutil

import fastaparser

########################################################################

if len(sys.argv) != 3:
	print 'FASTA-file converter from multi-line reads to one-line ones'
	print 'Usage: ', sys.argv[0], ' <input-file> <output-file>'
	sys.exit(0)

out = open(sys.argv[2], 'w')

fasta = fastaparser.read_fasta(sys.argv[1])
for name, seq in fasta:
    out.write(name + '\n')
    out.write(seq + '\n')

out.close()
