#!/usr/bin/env python
from __future__ import (print_function)

import argparse
import subprocess
import sys
import os
import os.path
import shutil
import yaml

from scripts.common import fill_default_values

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
parser.add_argument("--reuse-assemblies", type=str, help="Directory with existing assemblies to reuse")
parser.add_argument("--reuse-profiles", type=str, help="Directory with existing profiles to reuse")
parser.add_argument("--reuse-from", type=str, help="Directory with another assembly to reuse everything that is possible (overrides other --reuses)")
parser.add_argument("--no-stats", "-S", action="store_true", help="Skip the stats section (overrides the config value)")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")
parser.add_argument("--dryrun", action="store_true", help="Show tasks, do not execute them")

args = parser.parse_args()

exec_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
LOCAL_DIR = os.path.realpath(os.path.join(exec_dir, "../../../"))

base_params = ["snakemake", "--directory", os.path.realpath(args.dir), "--cores", str(args.threads), "--config", "LOCAL_DIR" + "=" + LOCAL_DIR]

if args.verbose:
    base_params.extend(["-p", "--verbose"])

if args.dryrun:
    base_params.extend(["--dryrun"])

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

print("Output folder set to", args.dir)

config_path = os.path.join(args.dir, "config.yaml")
if args.config:
    if os.path.exists(config_path):
        if subprocess.call(["diff", config_path, args.config]):
            print("\033[31mConfig path specified, but different config.yaml already exists in output folder", args.dir, "\033[0m")
            sys.exit(239)
    else:
        print("Copying config from", args.config)
        shutil.copy(args.config, config_path)

with cd(exec_dir):
    def call_snake(extra_params=[]):
        subprocess.check_call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)

    def reuse_dir(dir_from, dir_name):
        if not dir_from:
            return
        local_dir = os.path.join(args.dir, dir_name)
        if not os.path.isdir(dir_from):
            print("\033[33mWarning: {} source directory doesn't exist\033[0m".format(dir_from))
            return
        if os.path.exists(local_dir):
            print("\033[33mWarning: {} destination directory already exists\033[0m".format(dir_name))
            return
        os.symlink(dir_from, local_dir)

    with open(config_path) as config_in:
        config = yaml.load(config_in)
    fill_default_values(config)

    if args.reuse_from:
        args.reuse_assemblies = os.path.join(args.reuse_from, "assembly")
        args.reuse_profiles = os.path.join(args.reuse_from, "profile")

    reuse_dir(args.reuse_assemblies, "assembly")
    reuse_dir(args.reuse_profiles, "profile")

    print("Step #1 - Assembly")
    call_snake()

    if config.get("reassembly", dict()).get("enabled", True):
        print("Step #1b - Reassembly")
        call_snake(["--snakefile", "Reassembly.snake"])

    if not args.no_stats and len(config.get("stats", dict())) > 0:
        print("Step #2 - Stats")
        call_snake(["--snakefile", "Stats.snake"])
