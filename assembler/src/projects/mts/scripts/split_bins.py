#!/usr/bin/env python
from __future__ import print_function

import os
from os import path
import sys
from Bio import SeqIO
import common
import subprocess

def print_usage():
        print("Usage: split_bins.py <contigs> <binning info> <output directory> [-p]")

contigs = sys.argv[1]
sample, _ = path.splitext(path.basename(contigs))
out_dir = sys.argv[3]
prepend_name = False
if len(sys.argv) > 4 and sys.argv[4] == "-p":
    prepend_name = True

binning = common.load_annotation(sys.argv[2], False)

subprocess.call("rm -f {}/{}-*.fasta".format(out_dir, sample), shell=True)

cags = set()
for seq in SeqIO.parse(contigs, "fasta"):
    seq_id = seq.id
    if prepend_name:
        seq.id = sample + "-" + seq_id
        seq.description = ""
    for cag in binning.get(seq_id, []):
        filename = cag + ".fasta"
        if prepend_name:
            filename = sample + "-" + filename
        with open(path.join(out_dir, filename), "a") as output:
            SeqIO.write(seq, output, "fasta")
