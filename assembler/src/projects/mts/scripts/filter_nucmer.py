#!/usr/bin/python
from __future__ import print_function

import re
import sys

def print_usage():
    print("For a sample assembly aligned to a reference, outputs only contigs which were aligned more than <threshold> percent of their length total, and that percent.")
    print("Usage: filter_nucmer.py <nucmer_output> <threshold>")
    print("Parameters:")
    print("<nucmer_output> is typically contig_reports/nucmer_output/yoursample.coords.filtered from QUAST output")
    print("<threshold> is the minimal total alignment of a contig (0-100%)")

if len(sys.argv) != 3:
    print_usage()
    sys.exit(1)

with open(sys.argv[1]) as input:
    threshold = float(sys.argv[2])
    align_data = re.compile("\d+ \d+ \| \d+ \d+ \| \d+ (\d+) \| [\d.]+ \| [^ ]+ ([^\s]+length_(\d+)[^\s]+)")
    contig = None
    contig_len = 0
    align_len = 0
    def print_contig():
        per = 100.0 * align_len / contig_len
        if (per > threshold):
            print(contig, per)

    for line in input:
        res = align_data.search(line)
        if res is None:
            continue
        new_contig = res.group(2)
        if contig != new_contig:
            if contig is not None:
                print_contig()
            contig = new_contig
            contig_len = int(res.group(3))
            align_len = 0
        #Assuming that all alignments of the same contig are consequent
        align_len += int(res.group(1))
    #Print the last contig separately
    print_contig()
