#!/usr/bin/env python
from __future__ import (print_function)

from operator import itemgetter
#import subprocess
import os.path
import sys
import yaml

if len(sys.argv) < 3:
    print("Usage: choose_samples.py <input table> <input bins> <output table> <output dir> ")
    exit(1)

PROF = sys.argv[1]
BINS = sys.argv[2]
PROF_OUT = sys.argv[3]
DIR = sys.argv[4]
CAGS = set()
with open(BINS) as input:
    for line in input:
        bin = line.split("\t")[0]
        CAGS.add(bin)

DESIRED_ABUNDANCE = 999999 #sys.maxsize
MIN_ABUNDANCE = 4
MIN_TOTAL_ABUNDANCE = 15

prof_dict = dict()

make_excluded = True
excluded_dir = os.path.join(DIR, "excluded")

#Assuming that samples are enumerated consecutively from 1 to N
with open(PROF) as input:
    for line in input:
        params = line.split()
        CAG = params[0]
        if len(CAGS) == 0 and CAG not in CAGS:
            continue
        profile = list(map(float, params[1:]))

        print("Profile of", CAG, ":", profile)

        weighted_profile = list((i, ab)
            for i, ab in enumerate(profile) if ab >= MIN_ABUNDANCE) #and path.exists("{}/{}/sample{}_1.fastq".format(DIR, CAG, i + 1)))
        weighted_profile.sort(key = itemgetter(1))

        total = 0
        samples = []
        #If we have overabundant samples, use the least.
        try:
            i = next(x for x, _ in weighted_profile if profile[x] >= DESIRED_ABUNDANCE)
            total = profile[i]
            samples = [i + 1]
        except StopIteration:
            #If there isn't any, collect from samples, starting from the largest
            for i, _ in reversed(weighted_profile):
                total += profile[i]
                samples.append(i + 1)
                if total >= DESIRED_ABUNDANCE:
                    break

        print("Chosen samples are", samples, "with total mean abundance", sum)
        prof_dict[CAG] = total
        config_dir = DIR
        if total < MIN_TOTAL_ABUNDANCE:
            print(CAG, "is too scarce; skipping")
            if make_excluded and not os.path.isdir(excluded_dir):
                os.mkdir(excluded_dir)
            make_excluded = False
            config_dir = excluded_dir
        config_path = os.path.join(config_dir, CAG + ".info")
        with open(config_path, "w") as out:
            print("total", sum(profile), file=out)
            for i, ab in enumerate(profile, start=1):
                if i in samples:
                    print("+", file=out, end="")
                print("sample" + str(i), ab, file=out)

with open(PROF_OUT, "w") as prof_out:
    yaml.dump(prof_dict, prof_out)
