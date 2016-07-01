#!/usr/bin/python
from __future__ import print_function

import os
from os import path
import sys
from Bio import SeqIO
import common

def print_usage():
        print("Usage: split_bins.py <contigs> <binning info> <output directory>")

contigs = sys.argv[1]
sample, _, _ = path.basename(file).partition(".")
out_dir = sys.argv[3]

binning = common.load_annotation(sys.argv[2], False)

cags = set()
for seq in SeqIO.parse(contigs, "fasta"):
    #seq.id = common.get_id(seq.id, sample)
    #seq.description = ""
    for cag in binning.get(seq.id, []):
        with open(path.join(out_dir, "{}-{}.fasta".format(sample, cag)), "a") as output:
            SeqIO.write(seq, output, "fasta")
