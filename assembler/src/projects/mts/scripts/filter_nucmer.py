#!/usr/bin/python
from __future__ import print_function

import re
import sys
from os import path

def print_usage():
    print("For a sample assembly aligned to a reference, outputs only contigs which were aligned more than <threshold> percent of their length total, and that percent.")
    print("Usage: filter_nucmer.py <output> <length> <threshold>")
    print("Parameters:")
    print("<output> is the QUAST output directory. Files with contig names will be put there.")
    print("<length> is minimal contig length (default: INF)")
    print("<threshold> is the minimal total alignment of a contig (0-100%)")

if len(sys.argv) != 4:
    print_usage()
    sys.exit(1)

report_data = dict()
out_dir = sys.argv[1]
report = "{}/report.tsv".format(out_dir)
with open(report) as input:
    for line in input:
        vals = line.rstrip().split("\t")
        report_data[vals[0]] = vals[1:]

min_length = int(sys.argv[2])
threshold = float(sys.argv[3])

for (i, sample) in enumerate(report_data["Assembly"]):
    reference_length = int(report_data["Reference length"][i])
    nucmer_output = "{}/contigs_reports/nucmer_output/{}.coords.filtered".format(out_dir, sample)
    output = open("{}/{}.cont".format(out_dir, sample), "w")
    if not path.exists(nucmer_output):
        output.close()
        continue
    with open(nucmer_output) as input:
        align_data = re.compile("\d+ \d+ \| \d+ \d+ \| \d+ (\d+) \| [\d.]+ \| [^ ]+ NODE_(\d+)_length_(\d+)")
        contig = None
        total_len = 0
        contig_len = 0
        align_len = 0
        def process_contig():
            per = 100.0 * align_len / contig_len
            if per > threshold and contig_len >= min_length:
                print("{}-{}\t{}\t{}".format(sample, contig, contig_len, per), file=output)
                return align_len
            return 0
        for line in input:
            res = align_data.search(line)
            if res is None:
                continue
            new_contig = res.group(2)
            if contig != new_contig:
                if contig is not None:
                    total_len += process_contig()
                contig = new_contig
                contig_len = int(res.group(3))
                align_len = 0
            #Assuming that all alignments of the same contig are consequent
            align_len += int(res.group(1))
        #Print the last contig separately
        total_len += process_contig()
        print("total\t{}\t{}".format(total_len, 100.0 * total_len / reference_length), file=output)
    output.close()
