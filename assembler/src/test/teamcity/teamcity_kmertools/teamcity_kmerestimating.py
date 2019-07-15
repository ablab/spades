#!/usr/bin/python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing SPAdes
# provide a path to .yaml file with test description

import sys
import os
import shutil
import argparse
import subprocess
from traceback import print_exc
import filecmp

#Log class, use it, not print
class Log:

    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)
        sys.stdout.flush()

    def warn(self, s):
        msg = "WARNING: " + s + "\n"
        self.text += msg
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


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('info', metavar='CONFIG_FILE', type=str,  help='a path to .yaml file with test description')
    parser.add_argument("--path", "-p", help="custom directory to spades-kmer-estimating", type=str)
    args = parser.parse_args()
    return args


def load_info(info_filename):
    sys.path.append("./ext/src/python_libs/")
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml

    info = pyyaml.load(open(info_filename, 'r'))
    return info


# Create output folder
def create_output_dir(dataset_info):
    #make dirs and remembering history
    output_dir = os.path.join(dataset_info["output_dir"], dataset_info["name"])

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    return output_dir


# Save meta information about this teamcity.py run
def save_run_info(args, output_dir):
    run_info = open(os.path.join(output_dir, "test_run.info"), "w")
    run_info.write(".info file: " + args.info + "\n")
    if args.spades_path:
        run_info.write("path to spades-kmer-estimating: " + str(args.path) + "\n")
    run_info.close()


# Compile SPAdes
def compile_spades(args, dataset_info, working_dir):
    if 'spades_compile' not in dataset_info or dataset_info["spades_compile"]:
        comp_params = ' '
        if 'compilation_params' in dataset_info:
            comp_params = " ".join(dataset_info["compilation_params"])

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


def make_kmerestimating_cmd(args, dataset_info, working_dir, output_dir):
    cmd = os.path.join(working_dir, "spades-kmer-estimating")
    if "K" in dataset_info:
        cmd += " --kmer " + dataset_info["K"]

    cmd += " --dataset " + dataset_info["dataset_path"]
    cmd += " > " + os.path.join(output_dir, "log")
    return cmd

def parse_log(output_dir):
    with open(os.path.join(output_dir, "log")) as f:
        lines = f.readlines()
        kmernum = lines[-1].split()[-1]
    with open(os.path.join(output_dir, "log"), "w") as fw:
        fw.write(kmernum)

def run_kmerestimating(working_dir, args, dataset_info, output_dir):
    if args.path:
        working_dir = args.path
        log.log("Different spades-kmerestimating path specified: " + working_dir)

    cmd = make_kmerestimating_cmd(args, dataset_info, working_dir, output_dir)
    log.log("Launching: " + cmd)

    ecode = os.system(cmd)
    if ecode != 0:
        log.err("spades-kmer-estimating finished abnormally with exit code " + str(ecode))
        return 4
    parse_log(output_dir)
    return 0


def etalon_saves(dataset_info, output_dir):
    if 'etalon_saves' in dataset_info:
        log.log("Comparing etalon saves now")
        etalon_folder = dataset_info["etalon_saves"]

        dircmp = filecmp.dircmp(output_dir, os.path.join(etalon_folder))

        if dircmp.diff_files != []:
            log.err("Comparing etalon saves did not pass")
            return 12
    return 0


### main ###
try:
    if len(sys.argv) == 1:
        command = 'python {} -h'.format(sys.argv[0])
        subprocess.call(command, shell=True)
        sys.exit(1)

    sys.stderr = sys.stdout
    exit_code = 0
    args = parse_args()
    dataset_info = load_info(args.info)
    working_dir = os.getcwd()

    # compile
    ecode = compile_spades(args, dataset_info, working_dir)
    if ecode != 0:
        log.err("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(3)

    #run kmerestimating
    output_dir = create_output_dir(dataset_info)
    ecode = run_kmerestimating(os.path.join(working_dir, "build_spades", "bin"), args, dataset_info, output_dir)

    if (ecode != 0):
        sys.exit(ecode)

    #compare_etalon
    ecode = etalon_saves(dataset_info, output_dir)
    sys.exit(ecode)

except SystemExit:
    raise

except:
    log.err("The following unexpected error occured during the run:")
    print_exc()
    sys.exit(239)
