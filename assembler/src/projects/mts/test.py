#!/usr/bin/python
from __future__ import print_function

import argparse
import os
import os.path
import re
import shutil
import sys
import subprocess
from traceback import print_exc
import yaml

from scripts.common import Table

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s
        self.text += msg + "\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()

# Taken from teamcity.py
# Compile SPAdes
def compile_spades(args, dataset_info, working_dir):
    if not args.cfg_compilation:
        log.log("Forced to use current SPAdes build, will not compile SPAdes");
    elif 'spades_compile' not in dataset_info.__dict__ or dataset_info.spades_compile:
        comp_params = ' '
        if 'compilation_params' in dataset_info.__dict__:
            comp_params = " ".join(dataset_info.compilation_params)

        bin_dir = 'build_spades'
        if not os.path.exists(bin_dir):
            os.makedirs(bin_dir)
        os.chdir(bin_dir)

        #Compilation
        err_code = os.system('cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=' + working_dir + ' ' + os.path.join(working_dir, 'src') + comp_params)
        err_code = err_code | os.system('make -j 16')
        err_code = err_code | os.system('make install')

        os.chdir(working_dir)

        if err_code != 0:
            # Compile from the beginning if failed
            shutil.rmtree('bin', True)
            shutil.rmtree('build_spades', True)
            return os.system('./spades_compile.sh ' + comp_params)
    return 0

def compile_mts(workdir):
    #if not args.cfg_compilation:
    #    log.log("Forced to use current build, will not compile");
    #    return 0
    os.chdir(workdir)
    ecode = subprocess.call("./prepare_cfg")
    if ecode != 0:
        return ecode
    return subprocess.call(["make", "-C", "build/release/projects/mts"])

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="Config template")
    parser.add_argument("dir", help="Output directory")
    parser.add_argument("--saves", "-s", type=str)
    parser.add_argument("--no-clean", action="store_true")
    parser.add_argument("--etalons", "-e", type=str, help="Directory of GF etalons")
    args = parser.parse_args()
    return args

def prepare_config(args, workdir):
    with open(os.path.join(args.config)) as template:
        params = yaml.load(template)
        params["BIN"] = os.path.join(workdir, "build/release/bin")
        params["SCRIPTS"] = os.path.join(workdir, "src/projects/mts/scripts")
        with open(os.path.join(args.dir, "config.yaml"), "w") as config:
            config.write(yaml.dump(params))

def run_mts(args, workdir):
    if not args.no_clean:
        shutil.rmtree(args.dir, True)
    if not os.path.exists(args.dir):
        os.mkdir(args.dir)
        prepare_config(args, workdir)
    mts_args = ["./mts.py", "--stats", args.dir]
    if args.saves:
        log.log("Copying saves from" + args.saves)
        for saves_dir in ["assembly", "reassembly"]:
            full_dir = os.path.join(args.saves, saves_dir)
            if os.path.isdir(full_dir):
                #shutil.copytree(os.path.join(args.saves, saves_dir), os.path.join(args.dir, saves_dir))
                os.symlink(full_dir, os.path.join(args.dir, saves_dir))
            else:
                log.warn("No " + saves_dir + " dir provided; skipping")
        #Don't touch symlinked assemblies because it may corrupt other runs with the same dependencies
        #mts_args.append("--reuse-assemblies")
    os.chdir(os.path.join(workdir, "src/projects/mts"))
    return subprocess.call(mts_args)

def check_etalons(args, workdir):
    class mut:
        res = 0

    re_num = re.compile("-?\d+(?:\.\d+)?")
    def read_cell(str):
        maybe_num = re_num.search(str)
        if not maybe_num:
            return 0
        return float(maybe_num.group(0))

    #TODO: more configurable? Ideally, to set custom threshold for each cell

    #Margin values should stay close to margin, otherwise it's a false pos/neg
    pos_threshold = 95
    neg_threshold = 5
    #For the rest ("floating" clusters we're unsure of), we allow broader +/- margin
    threshold = 10

    def compare_gf(ref, cag, val1, val2):
        log.log("Comparing {} in {}: {} vs {}".format(cag, ref, val1, val2))
        et_val = read_cell(val1)
        est_val = read_cell(val2)
        lower = pos_threshold if et_val > pos_threshold else max(0,   et_val - threshold)
        upper = neg_threshold if et_val < neg_threshold else min(100, et_val + threshold)
        if est_val < lower:
            log.err("GF of {} in {} = {}% is less than expected {:.2f}%".format(cag, ref, est_val, lower))
            mut.res = 7
        elif est_val > upper:
            log.err("GF of {} in {} = {}% is higher than expected {:.2f}%".format(cag, ref, est_val, upper))
            mut.res = 7

    for file in os.listdir(args.etalons):
        etalon = os.path.join(args.etalons, file)
        estimated = os.path.join(args.dir, "stats", "summary", file)
        log.log("Trying to compare " + etalon + " and " + estimated)
        if not os.path.isfile(estimated):
            log.warn("No table provided for " + file)
            continue
        try:
            log.log("Loading " + etalon)
            et_table = Table.read(etalon, headers=True)
            log.log("Loading " + estimated)
            est_table = Table.read(estimated, headers=True)
            log.log("Comparing GF for " + file)
            et_table.zip_with(est_table, compare_gf)
        except:
            log.err("Cannot load {}".format(file))
            raise
    return mut.res

if __name__ == "__main__":
    try:
        sys.stderr = sys.stdout
        args = parse_args()
        workdir = os.getcwd()
        ecode = 0

        #compile
        #if compile_spades(args, dataset_info, working_dir) != 0:
        #    log.err("SPAdes compilation finished abnormally with exit code " + str(ecode))
        #    sys.exit(3)

        ecode = compile_mts(workdir)
        if ecode != 0:
            log.err("MTS compilation finished abnormally with exit code " + str(ecode))
            sys.exit(3)

        ecode = run_mts(args, workdir)
        if ecode != 0:
            log.err("Error while running MTS: " + str(ecode))

        if args.etalons:
            ecode = check_etalons(args, workdir)

        sys.exit(ecode)

    except SystemExit:
        raise

    except:
        log.err("The following unexpected error occured during the run:")
        print_exc()
        sys.exit(239)
