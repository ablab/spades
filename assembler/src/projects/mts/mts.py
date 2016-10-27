#!/usr/bin/python
from __future__ import (print_function)

import argparse
import subprocess
import sys
import os
import os.path

parser = argparse.ArgumentParser(description="MTS - Metagenomic Time Series")

parser.add_argument("--threads", "-t", type=int, default=8, help="Number of threads")
parser.add_argument("dir", type=str, help="Output directory where config.yaml is located")
parser.add_argument("--stats", "-s", action="store_true", help="Calculate stats (when the REFS parameter in config.yaml is provided)")
parser.add_argument("--reuse-assemblies", action="store_true", help="Use existing assemblies (put them in the corresponding folders)")

args = parser.parse_args()

LOCAL_DIR = os.path.realpath(os.path.join(os.getcwd(), "../../../"))

base_params = ["snakemake", "--directory", args.dir, "--cores", str(args.threads), "--config", "LOCAL_DIR" + "=" + LOCAL_DIR]

def call_snake(extra_params=[]):
    subprocess.check_call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)

print("Step #1 - Assembly")
if args.reuse_assemblies:
    call_snake(["assemble_all", "--touch"])
call_snake()

if args.stats:
    print("Step #2a - Assembly statistics")
    call_snake(["--snakefile", "Snakestats", "stats_all"])

    print("Step #2b - Reassembly statistics")
    call_snake(["--snakefile", "Snakestats", "stats_reassembly"])
