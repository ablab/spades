#!/usr/bin/python
from __future__ import print_function

import re
import sys
from os import path

def print_usage():
    print("For a sample assembly aligned to a reference, outputs only contigs which were aligned more than <threshold> percent of their length total, and that percent.")
    print("Usage: filter_nucmer.py <QUAST contig report> <length> <threshold>")
    print("Parameters:")
    print("<length> is minimal contig length")
    print("<threshold> is the minimal total alignment of a contig (0-100%)")

if len(sys.argv) != 4:
    print_usage()
    sys.exit(1)

nucmer_output_fn = sys.argv[1]
min_length = int(sys.argv[2])
threshold = float(sys.argv[3])

if not path.exists(nucmer_output_fn):
    print("File {} doesn't exist".format(nucmer_output_fn))
    sys.exit(2)

with open(nucmer_output_fn, "r") as nucmer_output:
    contig_data = re.compile("CONTIG: ([\w.-]+) \((\d+)bp\)")
    align_data = re.compile("Best alignment score: ([\d.]+)")
    split_format = re.compile("^([\w.-]+_)_(\d+_\d+)_$") #Replace brackets back
    while True:
        line = nucmer_output.readline()
        if not line or line.startswith("Analyzing coverage"):
            break
        res = contig_data.search(line)
        if res is None:
            continue
        contig = res.group(1)
        split = split_format.match(contig)
        if split:
            contig = split.group(1) + "(" + split.group(2) + ")"
        contig_len = int(res.group(2))
        line = nucmer_output.readline()
        if contig_len < min_length:
            continue #Too short
        res = align_data.search(line)
        if res is None:
            continue #Unaligned contig
        score = float(res.group(1))
        per = 100.0 * score / contig_len
        if per > threshold:
            print("{}\t{}\t{}".format(contig, contig_len, per))
