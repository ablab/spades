#!/usr/bin/env python
import os
import sys
import numpy as np
import pysam 
from Bio import SeqIO
import math

def binned_average(lens, bam_file, bin_size) :
    seq_bin_counts = [np.zeros(l / bin_size + 1) for l in lens]
    with pysam.AlignmentFile(bam_file, "rb") as sam:
        for read in sam:
            if not read.is_unmapped:
                seq_bin_counts[read.reference_id][read.reference_start / bin_size] += read.query_alignment_length
        
    return [a[:-1] / bin_size for a in seq_bin_counts]

if len(sys.argv) < 4: 
    print "Usage: %s <fragments fasta> <bam> <bin size>" % sys.argv[0]
    sys.exit(1) 

seq_infos = [(seq.name, len(seq.seq)) for seq in SeqIO.parse(open(sys.argv[1], "r"), "fasta")]

for (name, binned_avg) in zip([n for (n, l) in seq_infos], binned_average([l for (n, l) in seq_infos], sys.argv[2], int(sys.argv[3]))):
    for a in binned_avg:
        print name, a
