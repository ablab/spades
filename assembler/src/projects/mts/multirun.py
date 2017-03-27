#!/usr/bin/env python
from __future__ import (print_function)

import argparse
import os
import os.path
import subprocess

from scripts.common import load_dict, dump_dict

parser = argparse.ArgumentParser(description="MTS Multi Runner")

parser.add_argument("--threads", "-t", type=int, default=8, help="Number of threads for each run")
parser.add_argument("dir", type=str, help="Output directory")
parser.add_argument("--config", type=str, help="Base config")
#parser.add_argument("--in", type=str, help="Input directory")
#parser.add_argument("--ref", type=str, help="References directory")
parser.add_argument("--pipelines", type=str, nargs="+", default=["canopy"], help="Pipeline configurations to run")
parser.add_argument("--stats", "-s", action="store_true", help="Calculate stats (when the REFS parameter in config.yaml is provided)")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")

args = parser.parse_args()

with open(args.config) as config_in:
    config_template = load_dict(config_in)
assemblies_dir = dict()

for pipeline in args.pipelines:
    print("Running", pipeline)
    dir = os.path.join(args.dir, pipeline)
    if not os.path.exists(dir):
        os.makedirs(dir)
    call_params = ["./mts.py", "-t", str(args.threads), dir]
    if args.stats:
        call_params.extend(["--stats"])
    config = config_template.copy()
    params = pipeline.split("_")
    assembly_name = None
    if len(params) == 1: #Main pipeline: binner only
        config["BINNER"] = params[0]
        assembly_name = "main"
    else: #Alt pipeline: assembler + binner
        config["ASSEMBLER"] = params[0]
        config["BINNER"] = params[1]
        call_params.extend(["--alt"])
        assembly_name = params[0]
    with open(os.path.join(dir, "config.yaml"), "w") as config_out:
        dump_dict(config, config_out)
    # Try to reuse assemblies from previous runs with the same assembler
    assembly_dir = assemblies_dir.get(assembly_name)
    if assembly_dir:
        call_params.extend(["--reuse-assemblies", assembly_dir])
    else:
        assemblies_dir[assembly_name] = os.path.join(dir, "assembly")
    subprocess.check_call(call_params)

#TODO: compare stats
