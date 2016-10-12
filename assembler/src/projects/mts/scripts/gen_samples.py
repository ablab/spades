#!/usr/bin/env python3
from __future__ import print_function

import argparse
import os
import os.path
import random
import shutil
import subprocess
import sys
from common import gather_refs, dump_dict
from scipy.stats import expon

def gen_profile(args):
    if args.distribution == "uni":
        #def rand():
        #    return random.randint(0, args.scale)
        pass
    elif args.distribution == "exp":
        def rand():
            return int(expon.rvs(scale=args.scale))

    refs = dict(gather_refs(args.references))
    if args.dump_desc:
        with open(args.dump_desc, "w") as desc:
            dump_dict(refs, desc)
    for ref in refs:
        print(ref, end=" ")
        for _ in range(args.samples):
            print(rand(), end=" ")
        print()

def gen_samples(args):
    refs = dict(gather_refs(args.references.split(",")))
    try:
        os.mkdir(args.out_dir)
    except OSError:
        pass

    read_len = args.read_length
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
                ref_len = os.stat(ref_path).st_size
                reads = ref_len * abundance // read_len
                print("Generating", reads, "reads for subsample", i, "of", ref_name)
                sample_dir = os.path.join(args.out_dir, "sample" + str(i))
                if first_line:
                    shutil.rmtree(sample_dir, ignore_errors=True)
                    subprocess.check_call(["mkdir", "-p", sample_dir])

                temp_1 = sample_dir + ".tmp.r1.fastq"
                temp_2 = sample_dir + ".tmp.r2.fastq"
                subprocess.check_call(["wgsim", "-N", str(reads), "-r", "0", "-1", str(read_len), "-2", str(read_len), "-d", "300", "-s", "10", "-e", "{:.2f}".format(args.error_rate), "-S", str(i), ref_path, temp_1, temp_2], stdout=subprocess.DEVNULL)

                print("Merging temporary files")
                for temp, out in [(temp_1, os.path.join(sample_dir, "r1.fastq")), (temp_2, os.path.join(sample_dir, "r2.fastq"))]:
                    with open(temp) as input, open(out, "a") as output:
                        for line in input:
                            if line.startswith("IIIII"): #TODO: remove this hack
                                output.write(adj_qual)
                            else:
                                output.write(line)
                    os.remove(temp)
            print()
            first_line = False

parser = argparse.ArgumentParser(description="Metagenomic Time Series Simulator")
parser.add_argument("--references", "-r", type=str, help="Comma-separated list of references, or a directory with them, or a desc file with reference paths prepended with @", required=True)
subparsers = parser.add_subparsers()

gen_profile_args = subparsers.add_parser("prof", help="Generate a profile for the reference set")
gen_profile_args.add_argument("--dump-desc", "-d", type=str, help="Dump description file with reference paths")
gen_profile_args.add_argument("--samples", "-n", type=int, help="Sample count", default=1)
gen_profile_args.add_argument("--scale", "-s", type=int, help="Distribution scale", default=20)
gen_profile_args.add_argument("--distribution", "-t", choices=["uni", "exp"], help="Distribution type", default="uni")
gen_profile_args.set_defaults(func=gen_profile)

gen_samples_args = subparsers.add_parser("gen", help="Generate reads using a profile")
gen_samples_args.add_argument("--out-dir", "-o", type=str, help="Output directory. Will be totally overwritten!")
gen_samples_args.add_argument("--read-length", "-l", type=int, help="Read length", default=100)
gen_samples_args.add_argument("--error-rate", "-e", type=float, help="Base error rate", default=0)
gen_samples_args.add_argument("profile", type=str, help="File with reference profiles")
gen_samples_args.set_defaults(func=gen_samples)

args = parser.parse_args()
args.func(args)
