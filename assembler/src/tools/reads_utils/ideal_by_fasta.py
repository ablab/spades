#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import stat
import sys
import shutil

def read_fasta(filename):
    """
        Returns list of FASTA entries (in tuples: name, seq)
    """
    res_name = []
    res_seq = []
    first = True
    seq = ''

    for line in open(filename):
        if line[0] == '>':
            res_name.append(line.strip())
            if not first:
                res_seq.append(seq)
            else:
                first = False
            seq = ''
        else:
            seq += line.strip()
    res_seq.append(seq)
    return zip(res_name, res_seq)


def simulate_by_fasta(fasta, read_length, filename):
    outfile = open(filename, 'w')
    ideal = ''.join(['I' for i in range(0, read_length)])
    n = 0
    for name, seq in fasta:
        for i in xrange(0, len(seq) - read_length - 1):
            outfile.write('@' + name[1:] + '_' + str(n) + '_' + str(i) + '\n')
            outfile.write(seq[i : i + read_length] + '\n')
            outfile.write('+' + name[1:] + '_' + str(n) + '_' + str(i) + '\n')
            outfile.write(ideal + '\n')
        n += 1
    outfile.close()


simulate_by_fasta(read_fasta(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
