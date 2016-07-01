#!/usr/bin/python
from __future__ import print_function

import os
import shutil
import sys
import argparse
import subprocess
from traceback import print_exc
import yaml

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
    shutil.rmtree(args.dir, True)
    os.mkdir(args.dir)
    prepare_config(args, workdir)
    os.chdir(os.path.join(workdir, "src/projects/mts"))
    return subprocess.call(["./mts.py", "--stats", args.dir])

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

    sys.exit(ecode)

except:
    log.err("The following unexpected error occured during the run:")
    print_exc()
    sys.exit(239)
