#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import argparse
import glob
import logging
import os
import sys
from site import addsitedir
from os.path import abspath, dirname, realpath, join, isfile

python_modules_home = abspath(dirname(realpath(__file__)))
sys.path.append(join(python_modules_home, ".."))
import support


def remove_not_corrected_reads(output_dir):
    for not_corrected in glob.glob(os.path.join(output_dir, "*.bad.fastq")):
        os.remove(not_corrected)


def compress_dataset_files(input_file, ext_python_modules_home, max_threads, log, not_used_yaml_file, output_dir,
                           gzip_output):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith("2."):
        import pyyaml2 as pyyaml
        from joblib2 import Parallel, delayed
    elif sys.version.startswith("3."):
        import pyyaml3 as pyyaml
        from joblib3 import Parallel, delayed

    dataset_data = pyyaml.load(open(input_file))
    remove_not_corrected_reads(output_dir)
    is_changed = False
    if gzip_output:
        is_changed = True
        pigz_path = support.which("pigz")
        if pigz_path:
            compressor = "pigz"
        else:
            compressor = "gzip"
        log.info("\n== Compressing corrected reads (with %s)" % compressor)
        to_compress = []
        for reads_library in dataset_data:
            for key, value in reads_library.items():
                if key.endswith("reads"):
                    compressed_reads_filenames = []
                    for reads_file in value:
                        compressed_reads_filenames.append(reads_file + ".gz")
                        to_compress.append(reads_file)
                    reads_library[key] = compressed_reads_filenames

        log.info("\n== Files to compress: " + str(to_compress))
        if len(to_compress):
            for reads_file in to_compress:
                if not isfile(reads_file):
                    support.error(
                        "something went wrong and file with corrected reads (%s) is missing!" % reads_file, log)

            if pigz_path:
                for reads_file in to_compress:
                    support.sys_call([pigz_path, "-f", "-7", "-p", str(max_threads), reads_file], log)
            else:
                n_jobs = min(len(to_compress), max_threads)
                outputs = Parallel(n_jobs=n_jobs)(
                    delayed(support.sys_call)(["gzip", "-f", "-7", reads_file]) for reads_file in to_compress)
                for output in outputs:
                    if output:
                        log.info(output)
            log.info("\n== Files compression is finished")

    if not_used_yaml_file != "":
        is_changed = True
        not_used_dataset_data = pyyaml.load(open(not_used_yaml_file))
        dataset_data += not_used_dataset_data
        log.info("\n== Info about datasets not used in error correction stage is loaded")

    if is_changed:
        with open(input_file, 'w') as f:
            pyyaml.dump(dataset_data, f,
                        default_flow_style=False, default_style='"', width=float("inf"))
        log.info("\n== Dataset yaml file is updated")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",
                    help="path to input file",
                    action="store")
    parser.add_argument("--ext_python_modules_home",
                    help="path to ext python modules home",
                    action="store")
    parser.add_argument("--max_threads",
                    type=int,
                    help="max threads",
                    action="store")
    parser.add_argument("--output_dir",
                    help="path to output dir",
                    action="store")
    parser.add_argument("--gzip_output",
                    help="flag for enable gziping",
                    action="store_true")
    parser.add_argument("--not_used_yaml_file",
                    default="",
                    help="path to yaml file with not used data during error correction",
                    action="store")
    return parser.parse_args()


def main():
    args = parse_args()

    log = logging.getLogger("compressing")
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    compress_dataset_files(args.input_file, args.ext_python_modules_home, args.max_threads, log,
                           args.not_used_yaml_file, args.output_dir, args.gzip_output)


if __name__ == "__main__":
    main()
