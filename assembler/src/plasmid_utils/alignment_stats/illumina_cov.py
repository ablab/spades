#!/usr/bin/env python
from __future__ import print_function

import os
import sys
from pysam import AlignmentFile
from Bio import SeqIO, Seq
import math
read_len = 150
k = 99

def avg(a) :
    if len(a) == 0:
        return 0.
    return 1. * sum(a) / len(a)

def med(a) :
    if len(a) == 0:
        return 0.
    return a[len(a)/2]

def find_avg(a) :
    a.sort()
    step = len(a) / 10
    return (med(a))

def estimate(tslrs, sam_file, bin_size) :
    seqs = [len(rec.seq) for rec in SeqIO.parse(open(tslrs, "r"), "fasta")]
    sam = AlignmentFile(sam_file, "rb")
    covs = [[0.] * (seq_len / bin_size + 1) for seq_len in seqs]
    print('Processing %d long reads' % len(seqs))
    cnt = 0
    cur_limit = 0
    for read in sam:
        cnt += 1
        if cnt % 1000000 == 0:
            print(cnt)
        if not read.is_unmapped:
            covs[read.reference_id][read.reference_start / bin_size] += 1. * read.query_alignment_length / bin_size
    names = [[rec.name,len(str(rec.seq))] for rec in SeqIO.parse(open(tslrs, "r"), "fasta")]
    return zip([find_avg(a[:-1]) for a in covs], names)

if len(sys.argv) < 4: 
    print("Usage: %s <tslrs> <sam> <out_fn>" % sys.argv[0])
    sys.exit(1) 

tslrs = sys.argv[1]
sam_file = sys.argv[2]
with open(sys.argv[3], 'w') as out:
    out.write ("Name\tcov\tkmercov(check constants!)\tlength\n")
    for cov, name in sorted(estimate(tslrs, sam_file, 500), reverse=True):
        out.write("%s\t%f\t%f\t%d\n" % (name[0], cov, cov * (read_len-k)/ read_len, name[1]))
