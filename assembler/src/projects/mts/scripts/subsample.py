#!/usr/bin/env python
from __future__ import (print_function)

import argparse
import os
import os.path
import subprocess
import sys

parser = argparse.ArgumentParser(description="Subsample reads for highly-covered organisms")

parser.add_argument("file", type=str, help="File with chosen samples")
parser.add_argument("reads", type=str, help="Reads directory")
parser.add_argument("output", type=str, help="Output directory")

args = parser.parse_args()

MAX_COV = 100

with open(args.file) as samples_info:
    samples = []
    total = 0
    for line in samples_info:
        sample_data = line.split()
        if sample_data[0][0] == "+":
            sample = sample_data[0][1:]
            if not os.path.exists(os.path.join(args.reads, "{}_1.fastq.gz".format(sample))):
                print("\033[33mWarning: {} contains no reads\033[0m".format(sample))
                continue
            cov = float(sample_data[1])
            samples.append((sample, cov))
            total += cov

frac = MAX_COV / total
if frac > 1:
    frac = 1.0

subtotal = total * frac
print("Subsampling reads to reduce ", total, " to ", subtotal, "(", frac, "x)", sep="")
with open(os.path.join(args.output, "total.cov"), "w") as cov_out:
    print(subtotal, file=cov_out)

def run_subsampler(file, reads, frac):
    params = ["seqtk", "sample", "-s100", reads, str(frac)]
    print("Calling", " ".join(params))
    subprocess.check_call(params, stdout=file, stderr=sys.stdout)

left_file = open(os.path.join(args.output, "left.fq"), "w")
right_file = open(os.path.join(args.output, "right.fq"), "w")

for sample, cov in samples:
    for file, reads in [(left_file, os.path.join(args.reads, "{}_1.fastq.gz".format(sample))),
                        (right_file, os.path.join(args.reads, "{}_2.fastq.gz".format(sample)))]:
        run_subsampler(file, reads, frac)
