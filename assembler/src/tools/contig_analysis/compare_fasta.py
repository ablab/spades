#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys

def read_fasta(filename):
    """
        Returns set of FASTA strings
    """

    res_seq = []
    first = True
    seq = ''
    for line in open(filename):
        if line[0] == '>':
            if not first:
                res_seq.append(seq)
            else:
                first = False
            seq = ''
        else:
            seq += line.strip()
    res_seq.append(seq)
    return set(res_seq)


def compare_contigs(set1, set2):
	if len(set1) != len(set2):
		return False

	for contig in set1:
		if not contig in set2:
			return False

	return True


if len(sys.argv) != 3:
	print('Contigs comparator usage: ' + sys.argv[0] + ' <file1> <file2>')
	sys.exit(0)

set1 = read_fasta(sys.argv[1])
set2 = read_fasta(sys.argv[2])

print(compare_contigs(set1, set2))
