#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import argparse
import logging
import sys
import gzip
from os.path import abspath, dirname, realpath, join

# init python_lib_folder
python_modules_home = abspath(dirname(realpath(__file__)))
sys.path.append(join(python_modules_home, ".."))
import support


def write_single_read(in_file, out_file, read_name=None, is_fastq=False, is_python3=False):
    if read_name is None:
        read_name = support.process_readline(in_file.readline(), is_python3)
    if not read_name:
        return ""  # no next read
    read_value = support.process_readline(in_file.readline(), is_python3)
    line = support.process_readline(in_file.readline(), is_python3)
    fpos = in_file.tell()
    while (is_fastq and not line.startswith('+')) or (not is_fastq and not line.startswith('>')):
        read_value += line
        line = support.process_readline(in_file.readline(), is_python3)
        if not line:
            if fpos == in_file.tell():
                break
            fpos = in_file.tell()
    out_file.write(read_name + '\n')
    out_file.write(read_value + '\n')

    if is_fastq:
        read_quality = support.process_readline(in_file.readline(), is_python3)
        line = support.process_readline(in_file.readline(), is_python3)
        while not line.startswith('@'):
            read_quality += line
            line = support.process_readline(in_file.readline(), is_python3)
            if not line:
                if fpos == in_file.tell():
                    break
                fpos = in_file.tell()
        if len(read_value) != len(read_quality):
            support.error("The length of sequence and quality lines should be the same! "
                          "Check read %s (SEQ length is %d, QUAL length is %d)" %
                          (read_name, len(read_value), len(read_quality)))
        out_file.write("+\n")
        out_file.write(read_quality + '\n')
    return line  # next read name or empty string


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--args_filename",
                    help="path to file with args",
                    action="store")
    parser.add_argument("--dst",
                    help="path to dst dir",
                    action="store")
    return parser.parse_args()


def main():
    args = parse_args()

    # create logger
    log = logging.getLogger("Preprocess interlaced reads")
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    with open(args.args_filename) as f:
        lines = f.readlines()
        lines = [x.rstrip() for x in lines]
        for input_filename, out_left_filename, out_right_filename, was_compressed, is_fastq in \
                zip(lines[0::5], lines[1::5], lines[2::5], lines[3::5], lines[4::5]):
            was_compressed = (was_compressed == "True")
            is_fastq = (is_fastq == "True")

            if was_compressed:
                input_file = gzip.open(input_filename, 'r')
            else:
                input_file = open(input_filename)

            log.info("== Splitting %s into left and right reads (in %s directory)" % (input_filename, args.dst))
            out_files = [open(out_left_filename, 'w'), open(out_right_filename, 'w')]
            i = 0
            next_read_name = write_single_read(input_file, out_files[i], None, is_fastq,
                                               sys.version.startswith("3.") and was_compressed)
            while next_read_name:
                i = (i + 1) % 2
                next_read_name = write_single_read(input_file, out_files[i], next_read_name, is_fastq,
                                                   sys.version.startswith("3.") and was_compressed)
            if i == 0:
                support.error(
                    "the number of reads in file with interlaced reads (%s) should be EVEN!" % (input_filename),
                    log)
            out_files[0].close()
            out_files[1].close()

            input_file.close()


if __name__ == "__main__":
    main()
