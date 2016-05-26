#!/usr/bin/python
from __future__ import print_function

import os
from os import path
import sys
from Bio import SeqIO
import common

def print_usage():
        print("Usage: split_bins.py <sample> <binning info> <output directory>")

sample = sys.argv[1]
sample_name, _ = path.splitext(path.basename(sample))
out_dir = sys.argv[3]

binning = common.load_annotation(sample)

cags = set()
for seq in SeqIO.parse(sample, "fasta"):
    for cag in binning.get(common.get_id(seq.id), []):
        with open(path.join(out_dir, "{}-{}.fasta".format(sample_name, cag)), "a") as output:
            SeqIO.write(seq, output, "fasta")
