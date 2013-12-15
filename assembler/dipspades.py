#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import getopt

import spades
import support
import options_storage
import dipspades_logic

def command_line_to_str(argv):
    cl = ""
    for obj in argv:
        cl += str(obj) + " "
    return cl

def main():
    all_long_options = list(set(options_storage.long_options + dipspades_logic.DS_Args_List.long_options))
    all_short_options = options_storage.short_options + dipspades_logic.DS_Args_List.short_options

    dipspades_logic_py_command_line = "./dipspades_logic.py"
    spades_py_command_line = "./spades.py --diploid"

    try:
        options, not_options = getopt.gnu_getopt(sys.argv, all_short_options, all_long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        dipspades_logic.usage()
        sys.stderr.flush()
        sys.exit(1)
    if not options:
        dipspades_logic.usage()
        sys.stderr.flush()
        sys.exit(1)

    output_dir = None
    for opt, arg in options:
        # processing some special options
        if opt == '-o':
            #arg = os.path.abspath(arg)
            output_dir = os.path.abspath(arg) #arg
        elif opt == '--careful' or opt == '--mismatch-correction':
            continue
        if opt == '-h' or opt == '--help':
            dipspades_logic.usage()
            sys.stderr.flush()
            sys.exit(0)            
        # for all other options (and -o also)
        cur_opt_arg = " " + opt + " " + arg
        if opt.startswith("--"):  # long option
            if opt[2:] in options_storage.long_options or (opt[2:] + "=") in options_storage.long_options:
                spades_py_command_line += cur_opt_arg
            if opt[2:] in dipspades_logic.DS_Args_List.long_options or (opt[2:] + "=") in dipspades_logic.DS_Args_List.long_options:
                dipspades_logic_py_command_line += cur_opt_arg
        else: # short option
            if opt != '-o':
                if opt[1:] in options_storage.short_options:
                    spades_py_command_line += cur_opt_arg
                if opt[1:] in dipspades_logic.DS_Args_List.short_options:
                    dipspades_logic_py_command_line += cur_opt_arg

    if not output_dir:
        support.ds_error("The output_dir is not set! It is a mandatory parameter (-o output_dir).")

    spades_output_dir = os.path.join(output_dir, "spades")
    dipspades_output_dir = os.path.join(output_dir, "dipspades")

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(spades_output_dir):
        os.makedirs(spades_output_dir)
    if not os.path.isdir(dipspades_output_dir):
        os.makedirs(dipspades_output_dir)

    spades_py_command_line += " -o " + spades_output_dir
    dipspades_logic_py_command_line += " -o " + dipspades_output_dir
    #log = dipspades_logic.create_log(output_dir) TODO: log

    #print ("hello", options, not_options)
    #print (dipspades_logic_py_command_line)
    #print (spades_py_command_line)
    spades.main(spades_py_command_line.split())
    spades_result = os.path.join(spades_output_dir, "contigs.fasta")
    if not os.path.isfile(spades_result):
        support.ds_error("Something went wrong and SPAdes did not generate haplocontigs. DipSPAdes cannot proceed without them, aborting.")
    dipspades_logic_py_command_line += " --hap " + spades_result
    dipspades_logic.main(dipspades_logic_py_command_line, command_line_to_str(sys.argv), spades.spades_home, spades.bin_home)


if __name__ == '__main__':
    main()
