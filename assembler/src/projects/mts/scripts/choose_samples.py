#!/usr/bin/python
from __future__ import (print_function)

import glob
from operator import itemgetter
from os import path
import subprocess
import sys

if len(sys.argv) < 3:
    print("Usage: choose_samples.py <canopy.prof> <binning dir> [CAGS+]")
    exit(1)

PROF = sys.argv[1]
DIR = sys.argv[2]
CAGS = None
if len(sys.argv) == 4:
    CAGS = set(sys.argv[3:])
DESIRED_ABUNDANCE = 50
MIN_ABUNDANCE = 4
MIN_TOTAL_ABUNDANCE = 20

#Assuming that samples are enumerated consecutively from 1 to N
with open(PROF) as input:
    for line in input:
        params = line.split()
        CAG = params[0]
        if CAGS and CAG not in CAGS:
            continue
        profile = map(float, params[1:])

        print("Profile of", CAG, ":", profile)

        weighted_profile = list((i, ab)
            for i, ab in enumerate(profile) if ab >= MIN_ABUNDANCE and path.exists("{}/{}/sample{}_1.fastq".format(DIR, CAG, i + 1)))
        weighted_profile.sort(key = itemgetter(1))

        sum = 0
        samples = []
        #If we have overabundant samples, use the least.
        try:
            i = next(x for x, _ in weighted_profile if profile[x] >= DESIRED_ABUNDANCE)
            sum = profile[i]
            samples = [i + 1]
        except StopIteration:
            #If there isn't any, collect from samples, starting from the largest
            for i, _ in reversed(weighted_profile):
                sum += profile[i]
                samples.append(i + 1)
                if sum >= DESIRED_ABUNDANCE:
                    break

        print("Chosen samples are", samples, "with total mean abundance", sum)
        if sum < MIN_TOTAL_ABUNDANCE:
            print(CAG, "is too scarce; skipping")
            continue

        for suf, name in [("1", "left"), ("2", "right")]:
            reads = ["{}/{}/sample{}_{}.fastq".format(DIR, CAG, sample, suf) for sample in samples]
            with open("{}/{}/{}.fastq".format(DIR, CAG, name), "w") as output:
                subprocess.check_call(["cat"] + reads, stdout=output)
