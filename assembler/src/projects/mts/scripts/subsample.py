#!/usr/bin/env python
from __future__ import (print_function)

import argparse
import gzip
import os
import os.path
import subprocess
import sys

parser = argparse.ArgumentParser(description="Subsample reads for highly-covered organisms")

parser.add_argument("file", type=str, help=".info file with chosen samples")
parser.add_argument("reads", type=str, help="Input reads directory")
parser.add_argument("output", type=str, help="Output directory")

args = parser.parse_args()

MAX_COV = 100

def reads_file(sample, suffix=""):
    return os.path.join(args.reads, sample + suffix + ".fastq.gz")

with open(args.file) as samples_info:
    samples = []
    total = 0
    for line in samples_info:
        sample_data = line.split()
        if sample_data[0][0] == "+":
            sample = sample_data[0][1:]
            if not os.path.exists(reads_file(sample, "_1")):
                print("\033[33mWarning: {} contains no reads\033[0m".format(sample))
                continue
            cov = float(sample_data[1])
            samples.append(sample)
            total += cov

frac = MAX_COV / total
if frac > 1:
    frac = 1.0

output = open(os.path.join(args.output, "reads.info"), "w")

subtotal = total * frac
print(subtotal, file=output)

if frac < 1:
    print("Subsampling reads to reduce ", total, " to ", subtotal, "(", frac, "x)", sep="")

    filenames = tuple(os.path.join(args.output, basename + ".fq.gz") for name in ("left", "right"))

    def open_gzipped(filename):
        return subprocess.Popen(["gzip"], stdin=subprocess.PIPE, stdout=open(fullpath, "wb"))

    left_gz, right_gz = [open_gzipped(name) for name in filenames]

    def run_subsampler(file, reads, frac):
        params = ["seqtk", "sample", "-s100", reads, str(frac)]
        print("Calling", " ".join(params))
        subprocess.check_call(params, stdout=file.stdin, stderr=sys.stdout)

    for sample in samples:
        for file, reads in [(left_gz, os.path.join(args.reads, "{}_1.fastq.gz".format(sample))),
                            (right_gz, os.path.join(args.reads, "{}_2.fastq.gz".format(sample)))]:
            run_subsampler(file, reads, frac)

    files = [filenames]
else: #Nothing to subsample; reuse old files
    files = [(reads_file(sample, "_1"), reads_file(sample, "_2")) for sample in samples]

for left_file, right_file in files:
    print(left_file, right_file, file=output)
