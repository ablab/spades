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

binning = dict()
with open(sys.argv[2]) as binning_file:
    for line in binning_file:
        info = line.split(" : ")
        binning[common.get_id(info[0])] = info[1].split()
#         for cag in info[1].split():
#             if cag not in binning:
#                 binning[cag] = set()
#             binning[cag].add(info[0])

cags = set()
for seq in SeqIO.parse(sample, "fasta"):
    for cag in binning.get(common.get_id(seq.id), []):
        #if cag not in cags:
        #    os.mkdir(path.join(out_dir, cag))
        with open(path.join(out_dir, "{}-{}.fasta".format(sample_name, cag)), "a") as output:
            SeqIO.write(seq, output, "fasta")
