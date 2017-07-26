#!/usr/bin/env python
from __future__ import print_function
try:
    from itertools import izip as zip
except ImportError:
    pass

import re
import argparse
import os
import os.path
import sys

from common import contig_length

parser = argparse.ArgumentParser(description="Binner input formatter")

parser.add_argument("--type", "-t", choices=["binsanity", "canopy", "concoct", "gattaca", "maxbin"], help="Binner type", default="canopy")
parser.add_argument("--count", "-n", type=int, help="Number of data samples")
parser.add_argument("--jgi", action="store_true", help="Input data is in the JGI read-based profiler format")
parser.add_argument("profiles", type=str, help="Groups profiles in .tsv format")
parser.add_argument("out", type=str, help="Output path")

args = parser.parse_args()

class Formatter:
    def __init__(self, out, samples):
        self.out = open(out, "w")
        self.header(samples)

    def header(self, samples):
        pass

    def profile(self, contig, length, profile):
        print(contig, *profile, sep="\t", file=self.out)

class CanopyFormatter(Formatter):
    def profile(self, contig, length, profile):
        print(contig, *profile, file=self.out)

class ConcoctFormatter(Formatter):
    def header(self, samples):
        print(*(["contig"] + ["cov_mean_" + sample for sample in samples]), sep="\t", file=self.out)

class GattacaFormatter(Formatter):
    def header(self, samples):
        print(*(["contig", "length"] + ["cov_mean_" + sample for sample in samples]), sep="\t", file=self.out)

    def profile(self, contig, length, profile):
        print(contig, length, *profile, sep="\t", file=self.out)

class MaxBinFormatter:
    def __init__(self, out, samples):
        self.outs = list()
        with open(out, "w") as out_list:
            for sample in samples:
                name = out + "." + sample
                self.outs.append(open(name, "w"))
                print(name, file=out_list)

    def profile(self, contig, length, profile):
        for out, mpl in zip(self.outs, profile):
            print(contig, mpl, sep="\t", file=out)

formatters = {"binsanity": Formatter,
              "canopy": CanopyFormatter,
              "concoct": ConcoctFormatter,
              "gattaca": GattacaFormatter,
              "maxbin": MaxBinFormatter}

formatter = formatters[args.type](args.out, ["sample" + str(i) for i in range(1, args.count + 1)])

with open(args.profiles, "r") as input:
    if args.jgi:
        input.readline() #Skip the header from JGI format
    for line in input:
        params = line.strip().split("\t")
        length = params[1] if args.jgi else contig_length(params[0])
        values = params[3::2] if args.jgi else params[1:] #Skip variance columns from JGI format
        formatter.profile(params[0], length, values)
