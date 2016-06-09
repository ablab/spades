#!/usr/bin/python
from __future__ import (print_function)

import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(description="MTS - Metagenomic Time Series")

parser.add_argument("--threads", "-t", type=int, default=8, help="Number of threads")
parser.add_argument("dir", type=str, help="Output directory where config.yaml is located")
parser.add_argument("--stats", "-s", type=bool, default=False, help="Calculate stats (when the REFS parameter in config.yaml is provided)")

args = parser.parse_args()

base_params = ["snakemake", "--directory", args.dir, "--cores", str(args.threads)]
print(base_params)

def call_snake(extra_params=[]):
    subprocess.call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)

print("Step #1 - Assembly")
call_snake()

print("Step #2 - Bin filtering")
call_snake(["choose_all"])

print("Step #3 - Bin reassembly")
call_snake(["reassemble_all"])

if (args.stats):
    print("Step #4 - Statistics")
    call_snake(["stats_all"])
