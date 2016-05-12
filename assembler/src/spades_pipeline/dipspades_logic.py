#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import getopt
import os
import logging
import shutil
import errno
import options_storage
import support
import process_cfg
from distutils import dir_util
from os.path import abspath, expanduser


class DS_Args_List:
    long_options = "expect-gaps expect-rearrangements hap= threads= memory= tmp-dir= dsdebug hap-assembly dsK= saves= start-from=".split()
    short_options = "o:t:m:"


class DS_Args:
    max_threads = options_storage.THREADS
    max_memory = options_storage.MEMORY
    tmp_dir = None
    allow_gaps = False
    weak_align = False
    haplocontigs_fnames = []
    output_dir = ""
    haplocontigs = ""
    dev_mode = False
    haplotype_assembly = False
    k = 55
    saves = ""
    start_from = "dipspades"


def print_ds_args(ds_args, log):
    log.info("dipSPAdes parameters:")
    log.info("\tK value for dipSPAdes: " + str(ds_args.k))
    log.info("\tExpect gaps: " + str(ds_args.allow_gaps))
    log.info("\tExpect rearrangements: " + str(ds_args.weak_align))
    log.info("\tFiles with haplocontigs : " + str(ds_args.haplocontigs_fnames))
    log.info("\tHaplotype assembly stage: " + str(ds_args.haplotype_assembly))
    log.info("\tOutput directory: " + str(ds_args.output_dir))
    log.info("")
    log.info("\tDir for temp files: " + str(ds_args.tmp_dir))
    log.info("\tThreads: " + str(ds_args.max_threads))
    log.info("\tMemory limit (in Gb): " + str(ds_args.max_memory))


# src_config_dir - path of dipspades configs
def copy_configs(src_config_dir, dst_config_dir):
    if os.path.exists(dst_config_dir):
        shutil.rmtree(dst_config_dir)
    dir_util.copy_tree(src_config_dir, dst_config_dir, preserve_times=False)


def prepare_configs(src_config_dir, ds_args, log):
    config_dir = os.path.join(ds_args.output_dir, "dipspades_configs")
    copy_configs(src_config_dir, config_dir)
    #log.info("dipSPAdes configs were copied to " + config_dir)
    config_fname = os.path.join(config_dir, "config.info")
    return os.path.abspath(config_fname)


def write_haplocontigs_in_file(filename, haplocontigs):
    hapfile = open(filename, 'w')
    for hapcontig in haplocontigs:
        hapfile.write(hapcontig + "\n")
    hapfile.close()

def ParseStartPoint(start_point_arg, log):
    if start_point_arg == 'pbr':
        return 'dipspades:polymorphic_br'
    elif start_point_arg == 'kmg':
        return 'dipspades:kmer_gluer'
    elif start_point_arg == 'cc':
        return 'dipspades:consensus_construction'
    elif start_point_arg == 'ha':
        return 'dipspades:haplotype_assembly'
    log.info("ERROR: Start point " + start_point_arg + " was undefined")
    sys.exit(1)
        
def parse_arguments(argv, log):
    try:
        options, not_options = getopt.gnu_getopt(argv, DS_Args_List.short_options, DS_Args_List.long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        sys.stderr.flush()
        options_storage.usage("", dipspades=True)
        sys.exit(1)

    ds_args = DS_Args()
    for opt, arg in options:
        if opt == '-o':
            ds_args.output_dir = abspath(expanduser(arg))
        elif opt == '--expect-gaps':
            ds_args.allow_gaps = True
        elif opt == '--expect-rearrangements':
            ds_args.weak_align = True
        elif opt == '--hap':
            ds_args.haplocontigs_fnames.append(support.check_file_existence(arg, 'haplocontigs', log, dipspades=True))
        elif opt == '-t' or opt == "--threads":
            ds_args.max_threads = int(arg)
        elif opt == '-m' or opt == "--memory":
            ds_args.max_memory = int(arg)
        elif opt == '--tmp-dir':
            ds_args.tmp_dir = abspath(expanduser(arg))
        elif opt == '--dsdebug':
            ds_args.dev_mode = True
        elif opt == '--hap-assembly':
            ds_args.haplotype_assembly = True
        elif opt == '--dsK':
            ds_args.k = int(arg)
        elif opt == '--saves':
            ds_args.saves = os.path.abspath(arg)
            ds_args.dev_mode = True
        elif opt == '--start-from':
            ds_args.start_from = ParseStartPoint(arg, log)
            ds_args.dev_mode = True
    ds_args.haplocontigs = os.path.join(ds_args.output_dir, "haplocontigs")

    if not ds_args.output_dir:
        support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).", log, dipspades=True)
    if not ds_args.haplocontigs_fnames and ds_args.start_from == 'dipspades':
        support.error("cannot start dipSPAdes without at least one haplocontigs file!", log, dipspades=True)
    if not ds_args.tmp_dir:
        ds_args.tmp_dir = os.path.join(ds_args.output_dir, options_storage.TMP_DIR)

    if ds_args.start_from != 'dipspades' and ds_args.saves == '':
        support.error("saves were not defined! dipSPAdes can not start from " + ds_args.start_from)

    return ds_args


def prepare_config(config_fname, ds_args, log):
    args_dict = dict()
    args_dict["tails_lie_on_bulges"] = process_cfg.bool_to_str(not ds_args.allow_gaps)
    args_dict["align_bulge_sides"] = process_cfg.bool_to_str(not ds_args.weak_align)
    args_dict["haplocontigs"] = process_cfg.process_spaces(ds_args.haplocontigs)
    args_dict["output_dir"] = process_cfg.process_spaces(ds_args.output_dir)
    args_dict["developer_mode"] = process_cfg.bool_to_str(ds_args.dev_mode)
    args_dict["tmp_dir"] = process_cfg.process_spaces(ds_args.tmp_dir)
    args_dict["max_threads"] = ds_args.max_threads
    args_dict["max_memory"] = ds_args.max_memory
    args_dict["output_base"] = ""
    args_dict["ha_enabled"] = process_cfg.bool_to_str(ds_args.haplotype_assembly)
    args_dict["K"] = str(ds_args.k)
    args_dict['saves'] = ds_args.saves
    args_dict['entry_point'] = ds_args.start_from
    process_cfg.substitute_params(config_fname, args_dict, log)


def print_ds_output(output_dir, log):
    consensus_file = os.path.join(output_dir, "consensus_contigs.fasta")
    if os.path.exists(consensus_file):
        log.info(" * Assembled consensus contigs are in: " + support.process_spaces(consensus_file))

    paired_consensus_file = os.path.join(output_dir, "paired_consensus_contigs.fasta")
    if os.path.exists(paired_consensus_file):
        log.info(" * Assembled paired consensus contigs are in: " + support.process_spaces(paired_consensus_file))

    unpaired_consensus_file = os.path.join(output_dir, "unpaired_consensus_contigs.fasta")
    if os.path.exists(unpaired_consensus_file):
        log.info(" * Assembled unpaired consensus contigs are in: " + support.process_spaces(unpaired_consensus_file))

    hapalignment_file = os.path.join(output_dir, "haplocontigs_alignent")
    if os.path.exists(hapalignment_file):
        log.info(" * Alignment of haplocontigs is in: " + support.process_spaces(hapalignment_file))

    haplotype_assembly_file = os.path.join(output_dir, "haplotype_assembly.out")
    if os.path.exists(haplotype_assembly_file):
        log.info(" * Results of haplotype assembly are in: " + support.process_spaces(haplotype_assembly_file))

    consregions_file = os.path.join(output_dir, "conservative_regions.fasta")
    if os.path.exists(consregions_file):
        log.info(" * Conservative regions are in: " + support.process_spaces(consregions_file))

    possconsregions_file = os.path.join(output_dir, "possibly_conservative_regions.fasta")
    if os.path.exists(possconsregions_file):
        log.info(" * Possibly conservative regions are in: " + support.process_spaces(possconsregions_file))


def main(ds_args_list, general_args_list, spades_home, bin_home):
    log = logging.getLogger('dipspades')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    support.check_binaries(bin_home, log)
    ds_args = parse_arguments(ds_args_list, log)

    if not os.path.exists(ds_args.output_dir):
        os.makedirs(ds_args.output_dir)
    log_filename = os.path.join(ds_args.output_dir, "dipspades.log")
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)

    params_filename = os.path.join(ds_args.output_dir, "params.txt")
    params_handler = logging.FileHandler(params_filename, mode='a')
    log.addHandler(params_handler)

    log.info("\n")
    log.info("General command line: " + " ".join(general_args_list) + "\n")
    log.info("dipSPAdes command line: " + " ".join(ds_args_list) + "\n")
    print_ds_args(ds_args, log)
    log.removeHandler(params_handler)

    log.info("\n======= dipSPAdes started. Log can be found here: " + log_filename + "\n")
    write_haplocontigs_in_file(ds_args.haplocontigs, ds_args.haplocontigs_fnames)

    config_fname = prepare_configs(os.path.join(spades_home, "configs", "dipspades"), ds_args, log)
    ds_args.tmp_dir = support.get_tmp_dir(prefix="dipspades_", base_dir=ds_args.tmp_dir)
    prepare_config(config_fname, ds_args, log)

    try:
        log.info("===== Assembling started.\n")
        binary_path = os.path.join(bin_home, "dipspades")
        command = [binary_path, config_fname]
        support.sys_call(command, log)
        log.info("\n===== Assembling finished.\n")
        print_ds_output(ds_args.output_dir, log)
        if os.path.isdir(ds_args.tmp_dir):
            shutil.rmtree(ds_args.tmp_dir)
        log.info("\n======= dipSPAdes finished.\n")
        log.info("dipSPAdes log can be found here: " + log_filename + "\n")
        log.info("Thank you for using dipSPAdes!")
        log.removeHandler(log_handler)
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            if exc_type == OSError and exc_value.errno == errno.ENOEXEC: # Exec format error
                support.error("It looks like you are using SPAdes binaries for another platform.\n" +
                              support.get_spades_binaries_info_message(), dipspades=True)
            else:
                log.exception(exc_value)
                support.error("exception caught: %s" % exc_type, log)
    except BaseException: # since python 2.5 system-exiting exceptions (e.g. KeyboardInterrupt) are derived from BaseException
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            support.error("exception caught: %s" % exc_type, log, dipspades=True)


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
    main(sys.argv, "", spades_init.spades_home, spades_init.bin_home)
