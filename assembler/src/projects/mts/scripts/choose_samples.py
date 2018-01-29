#!/usr/bin/env python
from __future__ import (print_function)

import argparse
from operator import itemgetter
import os.path
import subprocess
import sys
import yaml

parser = argparse.ArgumentParser(description="Choose samples for bin reassembly")

parser.add_argument("-r", "--reads", type=str, help="Directory with reads")
parser.add_argument("-o", "--output", type=str, help="Output path")
parser.add_argument("prof", type=str, help="File with bin profiles")
parser.add_argument("-f", "--filtered", type=str, help="File with filtered bins")

args = parser.parse_args()

BINS = set()
if args.filtered:
    with open(args.filtered) as input:
        for line in input:
            bin = line.split()[0]
            BINS.add(bin)

MIN_TOTAL = 15
FAIR_TOTAL = 50
SOFT_MIN = 10
HARD_MIN = 4

prof_dict = dict()

make_excluded = True
excluded_dir = os.path.join(args.output, "excluded")

input = open(args.prof)
for line in input:
    exclude = False
    samples = []
    params = line.split()
    bin = params[0]
    profile = list(map(float, params[1:]))
    print("Profile of", bin, ":", profile)
    if BINS and bin not in BINS:
        print(bin, "was excluded from reassembly")
        exclude = True
    else:
        total = 0

        #Current strategy: choose all samples larger than soft threshold.
        #If there's not enough, use the lower hard threshold.
        #Sort samples by their abundancies
        weighted_profile = list((i, ab)
            for i, ab in enumerate(profile, start=1) if ab >= HARD_MIN)
        weighted_profile.sort(key=itemgetter(1), reverse=True)

        for i, cov in weighted_profile:
            if cov < SOFT_MIN and total >= FAIR_TOTAL: #Enough
                break
            if not os.path.exists(os.path.join(args.reads, bin, "sample{}_1.fastq.gz".format(i))):
                print("WARNING: sample", i, "does not contain reads for", bin)
                continue
            total += cov
            samples.append(i)

        print("Chosen samples are", samples, "with total mean abundance", total)
        prof_dict[bin] = total

        if total < MIN_TOTAL:
            print(bin, "is too scarce; excluding from default")
            exclude = True

    config_dir = args.output
    if exclude:
        if make_excluded and not os.path.isdir(excluded_dir):
            os.mkdir(excluded_dir)
        make_excluded = False
        config_dir = excluded_dir
    config_path = os.path.join(config_dir, bin + ".info")
    print("Dumping config into", config_path)
    with open(config_path, "w") as out:
        print("total", sum(profile), file=out)
        for i, ab in enumerate(profile, start=1):
            if i in samples:
                print("+", end="", file=out)
            print("sample" + str(i), ab, file=out)
