#!/usr/bin/env python
from __future__ import (print_function)

from operator import itemgetter
#import subprocess
import os.path
import sys
import yaml

if len(sys.argv) < 3:
    print("Usage: choose_samples.py <input table> <input bins> <input reads> <output dir> ")
    exit(1)

PROF = sys.argv[1]
FILTERED_BINS = sys.argv[2]
READ_DIR = sys.argv[3]
OUT_DIR = sys.argv[4]
BINS = set()
with open(FILTERED_BINS) as input:
    for line in input:
        bin = line.split()[0]
        BINS.add(bin)

DESIRED_ABUNDANCE = 50 #999999 #sys.maxsize
MIN_ABUNDANCE = 4
MIN_TOTAL_ABUNDANCE = 15

prof_dict = dict()

make_excluded = True
excluded_dir = os.path.join(OUT_DIR, "excluded")

#Assuming that samples are enumerated consecutively from 1 to N
#(it is forced by the pipeline)
with open(PROF) as input:
    for line in input:
        exclude = False
        samples = []
        params = line.split()
        bin = params[0]
        profile = list(map(float, params[1:]))
        print("Profile of", bin, ":", profile)
        if bin not in BINS:
            print(bin, "was excluded from reassembly")
            exclude = True
        else:
            #Sort samples by their abundancies
            weighted_profile = list((i, ab)
                for i, ab in enumerate(profile) if ab >= MIN_ABUNDANCE)
            weighted_profile.sort(key = itemgetter(1))

            total = 0
            #If we have overabundant samples, use the least.
            # try:
            #     i = next(x for x, _ in weighted_profile if profile[x] >= DESIRED_ABUNDANCE)
            #     total = profile[i]
            #     samples = [i + 1]
            # except StopIteration:

            #Current strategy: collect the desired abundance from samples, starting from the largest
            for i, _ in reversed(weighted_profile):
                if not os.path.exists(os.path.join(READ_DIR, bin, "sample{}_1.fastq.gz".format(i))):
                    print("WARNING: sample", i, "does not contain reads for", bin)
                    continue
                total += profile[i]
                samples.append(i + 1)
                if total >= DESIRED_ABUNDANCE:
                    break

            print("Chosen samples are", samples, "with total mean abundance", total)
            prof_dict[bin] = total

            if total < MIN_TOTAL_ABUNDANCE:
                print(bin, "is too scarce; skipping")
                exclude = True

        config_dir = OUT_DIR
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
