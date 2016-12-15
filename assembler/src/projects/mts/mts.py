#!/usr/bin/python
from __future__ import (print_function)

import argparse
import subprocess
import sys
import os
import os.path
import shutil

#copied from http://stackoverflow.com/questions/431684/how-do-i-cd-in-python/13197763#13197763
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

parser = argparse.ArgumentParser(description="MTS - Metagenomic Time Series")

parser.add_argument("--threads", "-t", type=int, default=8, help="Number of threads")
parser.add_argument("dir", type=str, help="Output directory")
parser.add_argument("--config", "-c", type=str, default="", help="config.yaml to be copied to the directory (unnecessary if config.yaml is already there)")
parser.add_argument("--stats", "-s", action="store_true", help="Calculate stats (when the REFS parameter in config.yaml is provided)")
parser.add_argument("--reuse-assemblies", action="store_true", help="Use existing assemblies (put them in the corresponding folders)")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")

args = parser.parse_args()

exec_dir=os.path.dirname(os.path.realpath(sys.argv[0]))
LOCAL_DIR = os.path.realpath(os.path.join(exec_dir, "../../../"))

base_params = ["snakemake", "--directory", os.path.realpath(args.dir), "--cores", str(args.threads), "--config", "LOCAL_DIR" + "=" + LOCAL_DIR]

if args.verbose:
    base_params.extend(["-p", "--verbose"])

if args.config:
    if os.path.exists(os.path.join(args.dir, "config.yaml")):
        print("Config path specified, but config.yaml already exists in output folder " + args.dir)
        sys.exit(239)

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

print("Output folder set to " + args.dir)

if args.config:
    print("Copying config from " + args.config)
    shutil.copy(args.config, args.dir)

with cd(exec_dir):
    def call_snake(extra_params=[]):
        subprocess.check_call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)
    
    print("Step #1 - Assembly")
    if args.reuse_assemblies:
        call_snake(["assemble_all", "--touch"])

    call_snake()
    
    if args.stats:
        print("Step #2a - Assembly statistics")
        call_snake(["--snakefile", "Stats.snake", "stats_all"])
    
        print("Step #2b - Reassembly statistics")
        call_snake(["--snakefile", "Stats.snake", "stats_reassembly"])

