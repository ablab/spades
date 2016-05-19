#!/usr/bin/python
from __future__ import print_function

import re
import sys

def print_usage():
    print("For a sample assembly aligned to a reference, outputs only contigs which were aligned more than <threshold> percent of their length total, and that percent.")
    print("Usage: filter_nucmer.py <quast_output> <sample> <threshold>")
    print("Parameters:")
    print("<quast_output> is the QUAST output directory")
    print("<sample> is the sample prefix (TODO: remove)")
    print("<threshold> is the minimal total alignment of a contig (0-100%)")

#if len(sys.argv) != 3:
if len(sys.argv) != 4:
    print_usage()
    sys.exit(1)

report = "{}/report.tsv".format(sys.argv[1])
with open(report) as input:
    report_data = dict(line.rstrip().split("\t", 1) for line in input)
sample = sys.argv[2]
#sample = report_data["Assembly"]
reference_length = int(report_data["Reference length"])
nucmer_output = "{}/contigs_reports/nucmer_output/{}.coords.filtered".format(sys.argv[1], report_data["Assembly"])
with open(nucmer_output) as input:
    #threshold = float(sys.argv[2])
    threshold = float(sys.argv[3])
    align_data = re.compile("\d+ \d+ \| \d+ \d+ \| \d+ (\d+) \| [\d.]+ \| [^ ]+ [^\s]+length_(\d+)[^\s]+ID_(\d+)")
    contig = None
    total_len = 0
    contig_len = 0
    align_len = 0
    def process_contig():
        per = 100.0 * align_len / contig_len
        if (per > threshold):
            print("{}-{}\t{}".format(sample, contig, per))
            return align_len
        return 0
    for line in input:
        res = align_data.search(line)
        if res is None:
            continue
        new_contig = res.group(3)
        if contig != new_contig:
            if contig is not None:
                total_len += process_contig()
            contig = new_contig
            contig_len = int(res.group(2))
            align_len = 0
        #Assuming that all alignments of the same contig are consequent
        align_len += int(res.group(1))
    #Print the last contig separately
    total_len += process_contig()
    print("total\t{}".format(100.0 * total_len / reference_length))
