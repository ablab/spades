#!/usr/bin/env python
from __future__ import print_function

from operator import itemgetter
import os
from os import path
import sys
from Bio import SeqIO
import common
import subprocess

def print_usage():
        print("Usage: split_bins.py <contigs> <binning info> <output directory> [-g]")

contigs = sys.argv[1]
sample, _ = path.splitext(path.basename(contigs))
out_dir = sys.argv[3]
glue = False
if len(sys.argv) > 4 and sys.argv[4]:
    glue = True

binning = common.load_annotation(sys.argv[2], False)

if glue:
    glue_binning = dict()
    for split, bins in binning.items():
        contig_bins = glue_binning.setdefault(common.extract_id(split), {})
        for bin in bins:
            contig_bins.setdefault(bin, 0)
            contig_bins[bin] += 1

if path.isdir(out_dir):
    subprocess.call("rm -f {}/*.fasta".format(out_dir), shell=True)
else:
    os.mkdir(out_dir)

for seq in SeqIO.parse(contigs, "fasta"):
    if glue:
        bins = []
        bin_freq = glue_binning.get(common.extract_id(seq.id))
        if bin_freq:
            bins = [max(bin_freq.items(), key=itemgetter(1))[0]]
    else:
        bins = binning.get(seq.id, [])
    for cag in bins:
        filename = cag + ".fasta"
        with open(path.join(out_dir, filename), "a") as output:
            SeqIO.write(seq, output, "fasta")
