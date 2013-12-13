#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import getopt
import os
import logging
import shutil

class DS_Args_List:
    long_options = "allow-gaps weak-align hap=".split()
    short_options = "o:"


class DS_Args:
    allow_gaps = False
    weak_align = False
    haplocontigs_fnames = []
    output_dir = ""
    haplocontigs = ""

    def print_fields(self):
        print ("allow_gaps: " + str(self.allow_gaps))
        print ("weak_align: " + str(self.weak_align))
        print ("files with haplocontigs : " + str(self.haplocontigs_fnames))
        print ("haplocontigs file: " + str(self.haplocontigs))
        print ("output_dir: " + str(self.output_dir))


def call_method(obj, method_name):
    getattr(obj, method_name)()


def usage():
    print ("./dipspades_logic.py: TODO: print usage")


def check_config_directory(output_base):
    config_dir = os.path.join(ds_args.output_dir, "configs")
    return config_dir

# src_config_dir - path of dipspades configs
def copy_configs(src_config_dir, dst_config_dir):
    if os.path.exists(dst_config_dir):
        shutil.rmtree(dst_config_dir)
    shutil.copytree(src_config_dir, dst_config_dir)


def prepare_configs(src_config_dir, ds_args, log):
    config_dir = os.path.join(ds_args.output_dir, "dipspades_configs")
    copy_configs(src_config_dir, config_dir)
    log.info("Configs were copied to " + config_dir)
    config_fname = os.path.join(config_dir, "config.info")
    return os.path.abspath(config_fname)


def write_haplocontigs_in_file(filename, haplocontigs):
    if os.path.exists(filename):
        os.remove(filename)
    hapfile = open(filename, 'a')
    for hapcontig in haplocontigs:
        hapfile.write(hapcontig + "\n")


def parse_arguments(argv):
    try:
        options, not_options = getopt.gnu_getopt(argv, DS_Args_List.short_options, DS_Args_List.long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        sys.stderr.flush()
        usage()
        sys.exit(1)

    ds_args = DS_Args()
    for opt, arg in options:
        if opt == '-o':
            ds_args.output_dir = arg
        elif opt == '--allow-gaps':
            ds_args.allow_gaps = True
        elif opt == '--weak-align':
            ds_args.weak_align = True
        elif opt == '--hap':
            ds_args.haplocontigs_fnames.append(arg)
    ds_args.haplocontigs = os.path.join(ds_args.output_dir, "haplocontigs")
    return ds_args


def check_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


def create_log(output_dir):
    log = logging.getLogger('dipspades')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    log_filename = os.path.join(output_dir, "dipspades.log")
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)
    return log


def check_binary(binary_dir, log):
    binary_path = os.path.join(binary_dir, "dipspades")
    if not os.path.isfile(binary_path):
        import support
        support.error("DipSPAdes binaries not found: " + binary_path, log)
    return binary_path


def get_dict_of_args(ds_args):
    import process_cfg
    args_dict = dict()
    args_dict["tails_lie_on_bulges"] = process_cfg.bool_to_str(ds_args.allow_gaps)
    args_dict["align_bulge_sides"] = process_cfg.bool_to_str(ds_args.weak_align)
    args_dict["haplocontigs"] = ds_args.haplocontigs
    args_dict["output_dir"] = ds_args.output_dir
    return args_dict


def prepare_config(config_fname, ds_args, log):
    import process_cfg
    args_dict = get_dict_of_args(ds_args)
    process_cfg.substitute_params(config_fname, args_dict, log)


def main(args, spades_home, bin_home):
    import support
    ds_args = parse_arguments(args)
    call_method(ds_args, "print_fields")

    check_output_dir(ds_args.output_dir)
    log = create_log(ds_args.output_dir)
    log.info("log was created")

    write_haplocontigs_in_file(ds_args.haplocontigs, ds_args.haplocontigs_fnames)
    log.info("file with haplocontigs was written\n")

    config_fname = prepare_configs(os.path.join(spades_home, "configs/dipspades"), ds_args, log)
    prepare_config(config_fname, ds_args, log)
    log.info("configs were prepared")

    binary_path = check_binary(bin_home, log)
    command = [binary_path, config_fname]
    support.sys_call(command, log)


if __name__ == '__main__':
    self_dir_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    spades_init_candidate1 = os.path.join(self_dir_path, "../../spades_init.py")
    spades_init_candidate2 = os.path.join(self_dir_path, "../../../bin/spades_init.py")
    if os.path.isfile(spades_init_candidate1):
        sys.path.append(os.path.dirname(spades_init_candidate1))
    elif os.path.isfile(spades_init_candidate2):
        sys.path.append(os.path.dirname(spades_init_candidate2))
    else:
        sys.stderr.write("Cannot find spades_init.py! Aborting..\n")
        sys.stderr.flush()
        sys.exit(1)
    import spades_init
    spades_init.init()
    main(sys.argv, spades_init.spades_home, spades_init.bin_home)
