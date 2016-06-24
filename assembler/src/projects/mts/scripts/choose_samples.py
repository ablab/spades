#!/usr/bin/python
from __future__ import (print_function)

import glob
from operator import itemgetter
from os import path
import subprocess
import sys

if len(sys.argv) != 4:
    print("Usage: choose_samples.py <CAG> <canopy.prof> <binning dir>")
    exit(1)

CAG = sys.argv[1]
PROF = sys.argv[2]
DIR = sys.argv[3]
DESIRED_ABUNDANCE = 50
MIN_ABUNDANCE = 4
MIN_TOTAL_ABUNDANCE = 20

#TODO: non-consecutive sample indexes
profile = [CAG] + map(float, subprocess.check_output(["grep", "-e", CAG + "\\s", PROF]).split()[1:])
print("Profile of", CAG, ":", profile)

weighted_profile = list((i, ab)
    for i, ab in enumerate(profile[1:], start=1) if ab >= MIN_ABUNDANCE and path.exists("{}/sample{}_1.fastq".format(DIR, i)))
weighted_profile.sort(key = itemgetter(1))
#print(weighted_profile)
sum = 0
samples = []
#If we have overabundant samples, use the least.
try:
    i= next(x for x, _ in weighted_profile if profile[x] >= DESIRED_ABUNDANCE)
    sum = profile[i]
    samples = [i]
except StopIteration:
    #If there isn't any, collect from samples, starting from the largest
    for i, _ in reversed(weighted_profile):
        sum += profile[i]
        samples.append(i)
        if sum >= DESIRED_ABUNDANCE:
            break

print("Chosen samples are", samples, "with total mean abundance", sum)
if sum < MIN_TOTAL_ABUNDANCE:
    print("The CAG is too scarce; skipping")
    exit(0)

for suf, name in [("1", "left"), ("2", "right")]:
    reads = ["{}/sample{}_{}.fastq".format(DIR, sample, suf) for sample in samples]
    with open("{}/{}.fastq".format(DIR, name), "w") as output:
        subprocess.check_call(["cat"] + reads, stdout=output)
