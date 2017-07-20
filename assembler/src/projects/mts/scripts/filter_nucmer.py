#!/usr/bin/python
from __future__ import print_function

import re
import sys
from os import path

def print_usage():
    print("For a sample assembly aligned to a reference, outputs only contigs which were aligned more than <threshold> percent of their length total, and that percent.")
    print("Usage: filter_nucmer.py <nucmer coords filtered> <length> <threshold>")
    print("Parameters:")
    print("<length> is minimal contig length (default: INF)")
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
    align_data = re.compile("\d+ \d+ \| \d+ \d+ \| \d+ (\d+) \| [\d.]+ \| [^ ]+ (NODE_(\d+)_length_(\d+)_.*$)")
    contig = None
    name = ""
    contig_len = 0
    align_len = 0
    def process_contig():
        per = 100.0 * align_len / contig_len
        if per > threshold and contig_len >= min_length:
            print("{}\t{}\t{}".format(name, contig_len, per))
            return align_len
        return 0
    for line in nucmer_output:
        res = align_data.search(line)
        if res is None:
            continue
        new_contig = res.group(3)
        if contig != new_contig:
            if contig is not None:
                process_contig()
            contig = new_contig
            name = res.group(2)
            contig_len = int(res.group(4))
            align_len = 0
        #Assuming that all alignments of the same contig are consequent
        align_len += int(res.group(1))
    #Print the last contig separately
    process_contig()
