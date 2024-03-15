#!/usr/bin/env python
from __future__ import (print_function)

import argparse
import os
import os.path
import subprocess
import sys
import yaml

parser = argparse.ArgumentParser(description="MTS Multi Runner")

all_assemblers = ["main", "spades", "megahit"]
all_binners = ["canopy", "concoct", "gattaca", "maxbin", "metabat"]
unsupported = set(["main_maxbin"])

parser.add_argument("--threads", "-t", type=int, default=8, help="Number of threads for each run")
parser.add_argument("dir", type=str, help="Output directory")
parser.add_argument("--config", "-c", type=str, default=None, help="Base config")
parser.add_argument("--pipelines", "-p", type=str, nargs="+", default=[], help="Pipeline configurations to run")
parser.add_argument("--assemblers", "-a", type=str, nargs="+", default=all_assemblers, help="Assemblers to use")
parser.add_argument("--binners", "-b", type=str, nargs="+", default=all_binners, help="Binners to use")
parser.add_argument("--exclude", "-e", type=str, nargs="+", default=[], help="Excluded (skipped) configurations")
parser.add_argument("--no-stats", "-S", action="store_true", help="Skip the stats section (overrides the config value)")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")
parser.add_argument("--reuse-from", type=str, default=None, help="Reuse assemblies from another multirun")
parser.add_argument("--ignore-errors", action="store_true")

args = parser.parse_args()

if not args.config:
    path = os.path.join(args.dir, "config.yaml")
    if not os.path.isfile(path):
        print("\033[31mError: no config provided with -c and no config found in the multirun directory\033[0m")
        sys.exit(1)
    args.config = path

with open(args.config) as config_in:
    config_template = yaml.load(config_in)

def pipelines():
    for assembler in args.assemblers:
        for binner in args.binners:
            yield assembler + "_" + binner
    for pipeline in args.pipelines:
        yield pipeline

prev_runs = dict()

excluded = unsupported.union(args.exclude)
for pipeline in pipelines():
    if pipeline in excluded:
        if pipeline in unsupported:
            print("\033[33mWarning:", pipeline, "is not currently supported; skipping\033[0m\n")
        continue
    print("Running", pipeline)
    cur_dir = os.path.join(args.dir, pipeline)
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)
    call_params = ["./mts.py", "-t", str(args.threads), cur_dir]
    if args.no_stats:
        call_params.extend(["--no-stats"])
    config = config_template.copy()
    for stage in ["assembly", "profile", "binning"]:
        config.setdefault(stage, dict())
    params = pipeline.split("_")
    assembly_name = params[0]
    if assembly_name == "main":
        config["profile"]["profiler"] = "mts"
    else:
        config["assembly"]["assembler"] = params[0]
        config["assembly"]["groups"] = ["*"]
        config["profile"]["profiler"] = "jgi"
        config["propagation"] = {"enabled": False}
        config["reassembly"] = {"enabled": False}

    config["binning"]["binner"] = params[1]
    with open(os.path.join(cur_dir, "config.yaml"), "w") as config_out:
        yaml.dump(config, config_out)
    # Try to reuse assemblies from previous runs with the same assembler
    prev_run = prev_runs.get(assembly_name)
    if prev_run:
        print("Reusing same data from", prev_run)
        call_params.extend(["--reuse-from", prev_run])
    elif args.reuse_from:
        for run in os.listdir(args.reuse_from):
            if run.startswith(assembly_name + "_"):
                path = os.path.join(args.reuse_from, run, "assembly")
                if os.path.isdir(path):
                    print("Reusing assembly from", path)
                    call_params.extend(["--reuse-assemblies", path])
                    break

    #TODO: rewrite using Snakemake API
    errcode = subprocess.call(call_params)
    if errcode:
        print(" ".join(call_params), "returned with error:", errcode)
        if not args.ignore_errors:
            sys.exit(errcode)
    elif not prev_run: #Reuse only successful run
        prev_runs[assembly_name] = cur_dir
    print()

#TODO: compare stats
