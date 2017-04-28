#!/usr/bin/env python
from __future__ import (print_function)

import argparse
import os
import os.path
import subprocess

from scripts.common import load_dict, dump_dict

parser = argparse.ArgumentParser(description="MTS Multi Runner")

all_assemblers = ["main", "spades", "megahit"]
all_binners = ["canopy", "concoct", "metabat"]
unsupported_configurations = set(["main_metabat", "spades_canopy", "megahit_canopy"])

parser.add_argument("--threads", "-t", type=int, default=8, help="Number of threads for each run")
parser.add_argument("dir", type=str, help="Output directory")
parser.add_argument("--config", "-c", type=str, help="Base config")
#parser.add_argument("--in", type=str, help="Input directory")
#parser.add_argument("--ref", type=str, help="References directory")
parser.add_argument("--pipelines", "-p", type=str, nargs="+", default=[], help="Pipeline configurations to run")
parser.add_argument("--assemblers", "-a", type=str, nargs="+", default=all_assemblers, help="Assemblers to use")
parser.add_argument("--binners", "-b", type=str, nargs="+", default=all_binners, help="Binners to use")
parser.add_argument("--exclude", "-e", type=str, nargs="+", default=[], help="Excluded (skipped) configurations")
parser.add_argument("--stats", "-s", action="store_true", help="Calculate stats (when the REFS parameter in config.yaml is provided)")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")

args = parser.parse_args()

with open(args.config) as config_in:
    config_template = load_dict(config_in)
assemblies_dir = dict()

def pipelines():
    for assembler in args.assemblers:
        for binner in args.binners:
            yield assembler + "_" + binner
    for pipeline in args.pipelines:
        yield pipeline

excluded = unsupported_configurations.union(args.exclude)
for pipeline in pipelines():
    if pipeline in excluded:
        if pipeline in unsupported_configurations:
            print("WARNING:", pipeline, "is not currently supported; skipping\n")
        continue
    print("Running", pipeline)
    dir = os.path.join(args.dir, pipeline)
    if not os.path.exists(dir):
        os.makedirs(dir)
    call_params = ["./mts.py", "-t", str(args.threads), dir]
    if args.stats:
        call_params.extend(["--stats"])
    config = config_template.copy()
    params = pipeline.split("_")
    assembly_name = params[0]
    if not assembly_name == "main":
        config["ASSEMBLER"] = params[0]
        call_params.extend(["--alt"])
    config["BINNER"] = params[1]
    with open(os.path.join(dir, "config.yaml"), "w") as config_out:
        dump_dict(config, config_out)
    # Try to reuse assemblies from previous runs with the same assembler
    assembly_dir = assemblies_dir.get(assembly_name)
    if assembly_dir:
        print("Reusing assemblies from", assembly_dir)
        call_params.extend(["--reuse-assemblies", assembly_dir])
    else:
        assemblies_dir[assembly_name] = os.path.join(dir, "assembly")
    errcode = subprocess.call(call_params)
    if errcode:
        print(" ".join(call_params), "returned with error:", errcode)
    print()

#TODO: compare stats
