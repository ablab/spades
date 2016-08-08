#!/usr/bin/python
from __future__ import (print_function)

import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(description="MTS - Metagenomic Time Series")

parser.add_argument("--threads", "-t", type=int, default=8, help="Number of threads")
parser.add_argument("dir", type=str, help="Output directory where config.yaml is located")
parser.add_argument("--stats", "-s", action="store_true", help="Calculate stats (when the REFS parameter in config.yaml is provided)")
parser.add_argument("--reuse-assemblies", action="store_true", help="Use existing assemblies (put them in the corresponding folders)")

args = parser.parse_args()

base_params = ["snakemake", "--directory", args.dir, "--cores", str(args.threads)]

def call_snake(extra_params=[]):
    subprocess.check_call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)

print("Step #1 - Assembly")
if args.reuse_assemblies:
    call_snake(["assemble_all", "--touch"])
call_snake()

print("Step #2 - Bin filtering")
call_snake(["choose_all"])

if args.stats:
    print("Step #2b - Assembly statistics")
    call_snake(["stats_all"])

print("Step #3 - Bin reassembly")
if args.reuse_assemblies and os.path.isdir(os.path.join(args.dir, "reassembly")):
    call_snake(["reassemble_all", "--reuse-assemblies"])
else:
    call_snake(["reassemble_all"])

if args.stats:
    print("Step #3b - Reassembly statistics")
    call_snake(["stats_reassembly"])
