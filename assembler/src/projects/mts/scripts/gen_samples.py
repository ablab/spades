#!/usr/bin/env python3
from __future__ import print_function

import argparse
import os
import os.path
import random
import shutil
import subprocess
import sys
import yaml
from common import gather_refs
#from scipy.stats import expon

def gen_profile(args):
    refs = dict(gather_refs(args.references))
    if args.dump_desc:
        with open(args.dump_desc, "w") as desc:
            yaml.dump(refs, desc)
    for ref in refs:
        print(ref, end=" ")
        for _ in range(args.samples):
            #abundance = expon.rvs(scale=30)
            abundance=random.randint(2, 20)
            print(int(abundance), end=" ")
        print()

def gen_samples(args):
    refs = dict(gather_refs(args.references))
    shutil.rmtree(args.out_dir, ignore_errors=True)
    os.mkdir(args.out_dir)

    read_len = 100
    adj_qual = "2" * read_len + "\n"

    with open(args.profile) as input:
        first_line = True
        for line in input:
            params = line.split()
            ref_name = params[0]
            ref_path = refs.get(ref_name)
            if not ref_path:
                print("Warning: no reference provided for", ref_name)
                continue
            for i, abundance in enumerate(map(int, params[1:]), start=1):
                ref_len = 5000000
                reads = ref_len * abundance // read_len
                print("Generating", reads, "reads for subsample", i, "of", ref_name)
                sample_dir = os.path.join(args.out_dir, "sample" + str(i))
                if first_line:
                    subprocess.check_call(["mkdir", "-p", sample_dir])

                temp_1 = sample_dir + ".tmp.r1.fastq"
                temp_2 = sample_dir + ".tmp.r2.fastq"
                subprocess.check_call(["wgsim", "-N", str(reads), "-r", "0.01", "-1", str(read_len), "-2", str(read_len), "-d", "3", "-s", "10", "-e", "0", "-S", str(i), ref_path, temp_1, temp_2], stdout=subprocess.DEVNULL)

                print("Merging temporary files")
                for temp, out in [(temp_1, os.path.join(sample_dir, "r1.fastq")), (temp_2, os.path.join(sample_dir, "r2.fastq"))]:
                    with open(temp) as input, open(out, "w") as output:
                        for line in input:
                            if line.startswith("IIIII"):
                                output.write(adj_qual)
                            else:
                                output.write(line)
                    os.remove(temp)
            print()
            first_line = False

parser = argparse.ArgumentParser(description="Metagenomic Time Series Simulator")
parser.add_argument("--references", "-r", type=str, help="List of references, or a directory with them, or a desc file with reference paths prepended with @", required=True)
subparsers = parser.add_subparsers()

gen_profile_args = subparsers.add_parser("prof", help="Generate a profile for the reference set")
gen_profile_args.add_argument("--dump-desc", "-d", type=str, help="Dump description file with reference paths")
gen_profile_args.add_argument("samples", type=int, help="Sample count")
gen_profile_args.set_defaults(func=gen_profile)

gen_samples_args = subparsers.add_parser("gen", help="Generate reads using a profile")
gen_samples_args.add_argument("--out-dir", "-o", type=str, help="Output directory", default="./")
gen_samples_args.add_argument("profile", type=str, help="File with reference profiles")
gen_samples_args.set_defaults(func=gen_samples)

args = parser.parse_args()
args.func(args)
