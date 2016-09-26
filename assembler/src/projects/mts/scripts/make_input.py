#!/usr/bin/env python
from __future__ import print_function
try:
    from itertools import izip as zip
except ImportError:
    pass

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Binner input formatter")
parser.add_argument("--type", "-t", type=str, help="Binner type (canopy or concoct)", default="canopy")
parser.add_argument("--output", "-o", type=str, help="Output file")
parser.add_argument("--dir", "-d", type=str, help="Directory with profiles (pairs of .id .mpl files)")
parser.add_argument("samples", type=str, nargs="+", help="Sample names")

args = parser.parse_args()

class CanopyFormatter:
    def __init__(self):
        pass

    def header(self, file, samples):
        pass

    def profile(self, file, contig, profile):
        print(contig, profile, file=out)

class ConcoctFormatter:
    def __init__(self):
        pass

    def header(self, file, samples):
        print("\t".join(["contig"] + ["cov_mean_" + sample for sample in samples]), file=out)

    def profile(self, file, contig, profile):
        print(contig.replace(",", "~"), profile.replace(" ", "\t"), sep="\t", file=out)

formatters = {"canopy": CanopyFormatter(), "concoct": ConcoctFormatter()}
formatter = formatters[args.type]

with open(args.output, "w") as out:
    formatter.header(out, args.samples)
    for sample in args.samples:
        id_file = "{}/{}.id".format(args.dir, sample)
        mpl_file = "{}/{}.mpl".format(args.dir, sample)

        print("Processing abundances from %s" % id_file)

        with open(id_file, "r") as ctg_id, open(mpl_file, "r") as ctg_mpl:
            for cid, cmpl in zip(ctg_id, ctg_mpl):
                formatter.profile(out, sample + "-" + cid.strip(), cmpl.strip())
