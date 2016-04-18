#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import getopt
from os.path import abspath, expanduser

import spades
import support
import options_storage
import dipspades_logic
import spades_init
spades_init.init()
spades_version = spades_init.spades_version


def main():
    all_long_options = list(set(options_storage.long_options + dipspades_logic.DS_Args_List.long_options))
    all_short_options = options_storage.short_options + dipspades_logic.DS_Args_List.short_options

    dipspades_logic_args = []
    spades_py_args = ["--diploid"]

    try:
        options, not_options = getopt.gnu_getopt(sys.argv, all_short_options, all_long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        options_storage.usage(spades_version, mode="dip")
        sys.stderr.flush()
        sys.exit(1)
    if not options:
        options_storage.usage(spades_version, mode="dip")
        sys.stderr.flush()
        sys.exit(1)

    output_dir = None
    spades_py_run_needed = False
    for opt, arg in options:
        # processing some special options
        if opt == '--test':
            output_dir = abspath("test_dipspades")
            spades_py_args = ["--diploid", "-1", os.path.join(spades_init.spades_home, "test_dataset/ecoli_1K_1.fq.gz"),
                              "-2", os.path.join(spades_init.spades_home, "test_dataset/ecoli_1K_2.fq.gz"), "--only-assembler"]
            dipspades_logic_args = []
            spades_py_run_needed = True
            break
        if opt == '-o':
            output_dir = abspath(expanduser(arg))
        elif opt == '--careful' or opt == '--mismatch-correction':
            continue
        if opt == '-v' or opt == '--version':
            options_storage.version(spades_version, mode="dip")
            sys.exit(0)
        if opt == '-h' or opt == '--help':
            options_storage.usage(spades_version, mode="dip")
            sys.exit(0)
        elif opt == "--help-hidden":
            options_storage.usage(spades_version, show_hidden=True, mode="dip")
            sys.exit(0)
        # for all other options
        cur_opt_arg = [opt]
        if arg:
            cur_opt_arg.append(arg)
        if opt.startswith("--"):  # long option
            if opt[2:] in options_storage.long_options or (opt[2:] + "=") in options_storage.long_options:
                spades_py_args += cur_opt_arg
                if opt[2:] in dipspades_logic.DS_Args_List.long_options or (opt[2:] + "=") in dipspades_logic.DS_Args_List.long_options:
                    dipspades_logic_args += cur_opt_arg
                else:
                    spades_py_run_needed = True
            else:
                dipspades_logic_args += cur_opt_arg
        else: # short option
            if opt != '-o':
                if opt[1:] in options_storage.short_options:
                    spades_py_args += cur_opt_arg
                    if opt[1:] in dipspades_logic.DS_Args_List.short_options:
                        dipspades_logic_args += cur_opt_arg
                    else:
                        spades_py_run_needed = True
                else:
                    dipspades_logic_args += cur_opt_arg

    if not output_dir:
        support.error("The output_dir is not set! It is a mandatory parameter (-o output_dir).", dipspades=True)

    spades_output_dir = os.path.join(output_dir, "spades")
    dipspades_output_dir = os.path.join(output_dir, "dipspades")

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(spades_output_dir):
        os.makedirs(spades_output_dir)
    if not os.path.isdir(dipspades_output_dir):
        os.makedirs(dipspades_output_dir)

    spades_result = ""
    if spades_py_run_needed:
        spades_py_args += ["-o", spades_output_dir]
        spades.main(spades_py_args)
        spades_result = os.path.join(spades_output_dir, "contigs.fasta")
        if not os.path.isfile(spades_result):
            support.error("Something went wrong and SPAdes did not generate haplocontigs. "
                      "DipSPAdes cannot proceed without them, aborting.", dipspades=True)

    dipspades_logic_args += ["-o", dipspades_output_dir]
    if spades_result != "":
        dipspades_logic_args += ["--hap", spades_result]
    dipspades_logic.main(dipspades_logic_args, sys.argv, spades.spades_home, spades.bin_home)


if __name__ == '__main__':
    main()
