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

parser.add_argument("--type", "-t", choices=["canopy", "concoct", "gattaca", "binsanity"], help="Binner type", default="canopy")
parser.add_argument("--count", "-n", type=int, help="Number of data samples")
parser.add_argument("--output", "-o", type=str, help="Output file")
parser.add_argument("profiles", type=str, help="Groups profiles in .tsv format")

args = parser.parse_args()

class CanopyFormatter:
    def __init__(self):
        pass

    def header(self, file, samples):
        pass

    def profile(self, file, contig, profile):
        print(contig, " ".join(profile), file=out)

class ConcoctFormatter:
    def __init__(self):
        pass

    def header(self, file, samples):
        print("\t".join(["contig"] + ["cov_mean_" + sample for sample in samples]), file=out)

    def profile(self, file, contig, profile):
        print(contig, *profile, sep="\t", file=out)

class BinSanityFormatter:
    def __init__(self):
        pass

    def header(self, file, samples):
        pass

    def profile(self, file, contig, profile):
        print(contig, *profile, sep="\t", file=out)

class GattacaFormatter:
    def __init__(self):
        pass

    def header(self, file, samples):
        print("\t".join(["contig", "length"] + ["cov_mean_" + sample for sample in samples]), file=out)

    def profile(self, file, contig, profile):
        l = contig_length(contig)
        print(contig, l, *profile, sep="\t", file=out)

formatters = {"canopy": CanopyFormatter(), "concoct": ConcoctFormatter(), "gattaca": GattacaFormatter(), "binsanity": BinSanityFormatter()}
formatter = formatters[args.type]

with open(args.output, "w") as out:
    formatter.header(out, ["sample" + str(i) for i in range(1, args.count + 1)])
    with open(args.profiles, "r") as input:
        for line in input:
            params = line.strip().split("\t")
            formatter.profile(out, params[0], params[1:])
