#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
from Bio import SeqIO

#### for MODE 1 ####
# deletion
del_start = 20000
del_len = 1500
# copy
copy_start = 40000
copy_len = 500
# relocation
rel_start = 60000
rel_len = 2000
rel_new_start = 80000
####

#### for MODE 2 ####
trans_to = 50000
trans_from = 40000
trans_len = 10000
####

if len(sys.argv) < 2:
    print('Usage: ' + sys.argv[0] + ' reference.fasta [reference2.fasta]')

in_fpath = sys.argv[1]
original = SeqIO.read(open(in_fpath), "fasta")
modified = original

if len(sys.argv) == 2:  # MODE1
    out_fpath = os.path.splitext(in_fpath)[0] + '.sv' + os.path.splitext(in_fpath)[1]
    modified._seq = modified.seq[:rel_start] + modified.seq[rel_start + rel_len:rel_new_start] + modified.seq[rel_start:rel_start + rel_len] + modified.seq[rel_new_start]
    modified._seq = modified.seq[:copy_start + copy_len] + modified.seq[copy_start:copy_start + copy_len] + modified.seq[copy_start + copy_len:]
    modified._seq = modified.seq[:del_start] + modified.seq[del_start + del_len:]

else:  # MODE2
    out_fpath = os.path.splitext(in_fpath)[0] + '.trans' + os.path.splitext(in_fpath)[1]
    in2_fpath = sys.argv[2]
    original2 = SeqIO.read(open(in2_fpath), "fasta")
    modified._seq = original.seq[:trans_to] + original2.seq[trans_from: trans_from + trans_len] +  original.seq[trans_to:]

output_handle = open(out_fpath, "w")
SeqIO.write([modified], output_handle, "fasta")
output_handle.close()

