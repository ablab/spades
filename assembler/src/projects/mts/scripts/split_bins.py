#!/usr/bin/python
from __future__ import print_function

import os
from os import path
import sys
from Bio import SeqIO
import common
import subprocess

def print_usage():
        print("Usage: split_bins.py <contigs> <binning info> <output directory>")

contigs = sys.argv[1]
sample, _ = path.splitext(path.basename(contigs))
out_dir = sys.argv[3]

binning = common.load_annotation(sys.argv[2], False)

subprocess.call("rm -f {}/{}-*.fasta".format(out_dir, sample), shell=True)

cags = set()
for seq in SeqIO.parse(contigs, "fasta"):
    seq_id = seq.id
    seq.id = sample + "-" + seq_id
    #seq.id = common.get_id(seq.id, sample)
    seq.description = ""
    for cag in binning.get(seq_id, []):
        with open(path.join(out_dir, "{}-{}.fasta".format(sample, cag)), "a") as output:
            SeqIO.write(seq, output, "fasta")
