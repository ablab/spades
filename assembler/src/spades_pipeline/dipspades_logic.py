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
import options_storage

class DS_Args_List:
    long_options = "expect-gaps expect-rearrangements hap=".split()
    short_options = "o:"


class DS_Args:
    allow_gaps = False
    weak_align = False
    haplocontigs_fnames = []
    output_dir = ""
    haplocontigs = ""

def print_ds_args(ds_args, log):
    log.info("dipSPAdes parameters:")
    log.info("\tExpect gaps: " + str(ds_args.allow_gaps))
    log.info("\tExpect rearrengements: " + str(ds_args.weak_align))
    log.info("\tFiles with haplocontigs : " + str(ds_args.haplocontigs_fnames))
    log.info("\tOutput directory: " + str(ds_args.output_dir))


def call_method(obj, method_name):
    getattr(obj, method_name)()


def usage():
    sys.stderr.write("dipSPAdes 1.0: genome assembler designed for diploid genomes with high heterozygosity rate\n")
    sys.stderr.write("Usage: " + str(sys.argv[0]) + " [options] -o <output_dir>" + "\n")
    sys.stderr.write("" + "\n")
    sys.stderr.write("Basic options:" + "\n")
    sys.stderr.write("-o\t<output_dir>\tdirectory to store all the resulting files (required)" + "\n")
    sys.stderr.write("--test\t\t\truns SPAdes on toy dataset" + "\n")
    sys.stderr.write("-h/--help\t\tprints this usage message" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Input reads:" + "\n")
    sys.stderr.write("--12\t<filename>\tfile with interlaced forward and reverse"\
                         " paired-end reads" + "\n")
    sys.stderr.write("-1\t<filename>\tfile with forward paired-end reads" + "\n")
    sys.stderr.write("-2\t<filename>\tfile with reverse paired-end reads" + "\n")
    sys.stderr.write("-s\t<filename>\tfile with unpaired reads" + "\n")
    sys.stderr.write("--pe<#>-12\t<filename>\tfile with interlaced"\
                         " reads for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-1\t<filename>\tfile with forward reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-2\t<filename>\tfile with reverse reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-s\t<filename>\tfile with unpaired reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-<or>\torientation of reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)" + "\n")
    sys.stderr.write("--mp<#>-12\t<filename>\tfile with interlaced"\
                         " reads for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-1\t<filename>\tfile with forward reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-2\t<filename>\tfile with reverse reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-s\t<filename>\tfile with unpaired reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-<or>\torientation of reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)" + "\n")
    sys.stderr.write("--sanger\t<filename>\tfile with Sanger reads\n")
    sys.stderr.write("--pacbio\t<filename>\tfile with PacBio reads\n")
    sys.stderr.write("--trusted-contigs\t<filename>\tfile with trusted contigs\n")
    sys.stderr.write("--untrusted-contigs\t<filename>\tfile with untrusted contigs\n")
    sys.stderr.write("Input haplocontigs:" + "\n")
    sys.stderr.write("--hap\t<filename>\tfile with haplocontigs" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Pipeline options:" + "\n")
    sys.stderr.write("--only-assembler\truns only assembling (without read error"\
                         " correction)" + "\n")
    sys.stderr.write("--disable-gzip-output\tforces error correction not to"\
                         " compress the corrected reads" + "\n")
    sys.stderr.write("--disable-rr\t\tdisables repeat resolution stage"\
                     " of assembling" + "\n")
                     
    sys.stderr.write("" + "\n")
    sys.stderr.write("DipSPAdes options:" + "\n")
    sys.stderr.write("--expect-gaps\tindicate that significant number of gaps in coverage is expected" + "\n")
    sys.stderr.write("--expect-rearrangements\tindicate that significant number of rearrangngements between haplomes of diploid genome is expected" + "\n")
    
    sys.stderr.write("" + "\n")
    sys.stderr.write("Advanced options:" + "\n")
    sys.stderr.write("--dataset\t<filename>\tfile with dataset description in YAML format" + "\n")
    sys.stderr.write("-t/--threads\t<int>\t\tnumber of threads" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % options_storage.THREADS)
    sys.stderr.write("-m/--memory\t<int>\t\tRAM limit for SPAdes in Gb"\
                         " (terminates if exceeded)" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % options_storage.MEMORY)
    sys.stderr.write("--tmp-dir\t<dirname>\tdirectory for temporary files" + "\n")
    sys.stderr.write("\t\t\t\t[default: <output_dir>/tmp]" + "\n")
    sys.stderr.write("-k\t\t<int,int,...>\tcomma-separated list of k-mer sizes"\
                         " (must be odd and" + "\n")
    sys.stderr.write("\t\t\t\tless than " + str(options_storage.MAX_K + 1) + ") [default: 'auto']" + "\n") # ",".join(map(str, k_mers_short))
    sys.stderr.write("--phred-offset\t<33 or 64>\tPHRED quality offset in the"\
                         " input reads (33 or 64)" + "\n")
    sys.stderr.write("\t\t\t\t[default: auto-detect]" + "\n")
    sys.stderr.write("--debug\t\t\t\truns SPAdes in debug mode (keeps intermediate output)" + "\n")


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
    #log.info("dipSPAdes configs were copied to " + config_dir)
    config_fname = os.path.join(config_dir, "config.info")
    return os.path.abspath(config_fname)


def write_haplocontigs_in_file(filename, haplocontigs):
    if os.path.exists(filename):
        os.remove(filename)
    hapfile = open(filename, 'a')
    for hapcontig in haplocontigs:
        if not os.path.exists(hapcontig):
            sys.stderr.write(hapcontig + ": file is nor found\n")
            sys.stderr.flush()
            sys.exit(1)
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
        elif opt == '--expect-gaps':
            ds_args.allow_gaps = True
        elif opt == '--expect-rearrangements':
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
        support.ds_error("DipSPAdes binaries not found: " + binary_path, log)
    return binary_path


def get_dict_of_args(ds_args):
    import process_cfg
    args_dict = dict()
    args_dict["tails_lie_on_bulges"] = process_cfg.bool_to_str(ds_args.allow_gaps)
    args_dict["align_bulge_sides"] = process_cfg.bool_to_str(ds_args.weak_align)
    args_dict["haplocontigs"] = ds_args.haplocontigs
    args_dict["output_dir"] = ds_args.output_dir
    args_dict["developer_mode"] = "false" #process_cfg.bool_to_str(False)
    return args_dict


def prepare_config(config_fname, ds_args, log):
    import process_cfg
    args_dict = get_dict_of_args(ds_args)
    process_cfg.substitute_params(config_fname, args_dict, log)

def print_ds_output(output_dir, log):
    consensus_file = os.path.join(output_dir, "consensus_contigs.fasta")
    if os.path.exists(consensus_file):
        log.info(" * Assembled consensus contigs are in: " + consensus_file)

    paired_consensus_file = os.path.join(output_dir, "paired_consensus_contigs.fasta")
    if os.path.exists(paired_consensus_file):
        log.info(" * Assembled paired consensus contigs are in: " + paired_consensus_file)

    unpaired_consensus_file = os.path.join(output_dir, "unpaired_consensus_contigs.fasta")
    if os.path.exists(unpaired_consensus_file):
        log.info(" * Assembled paired consensus contigs are in: " + unpaired_consensus_file)

    hapalignment_file = os.path.join(output_dir, "haplocontigs_alignent")
    if os.path.exists(hapalignment_file):
        log.info(" * Alignment of haplocontigs is in: " + hapalignment_file)

    haplotype_assembly_file = os.path.join(output_dir, "haplotype_assembly.out")
    if os.path.exists(haplotype_assembly_file):
        log.info(" * Assembled paired consensus contigs are in: " + haplotype_assembly_file)

    consregions_file = os.path.join(output_dir, "conservative_regions.fasta")
    if os.path.exists(consregions_file):
        log.info(" * Conservative regions are in: " + consregions_file)

    possconsregions_file = os.path.join(output_dir, "possibly_conservative_regions.fasta")
    if os.path.exists(possconsregions_file):
        log.info(" * Possibly conservative regions are in: " + possconsregions_file)

def write_params(output_dir, command_line, ds_command_line):
    params = os.path.join(output_dir, "params.txt")
    if os.path.exists(params):
        os.remove(params)
    params_file = open(params, 'a')
    params_file.write("Command line: " + command_line + "\n")
    params_file.write("dipSPAdes command line: " + ds_command_line + "\n")

def main(ds_command_line, general_command_line, spades_home, bin_home):
    import support
    
    args = ds_command_line.split()
    ds_args = parse_arguments(args)

    check_output_dir(ds_args.output_dir)
    log = create_log(ds_args.output_dir)

    log.info("\n")
    log.info("Command line: " + general_command_line + "\n")
    log.info("dipSPAdes command line: "+ ds_command_line + "\n")
    print_ds_args(ds_args, log)
    log.info("\n======= dipSPAdes started. Log can be found here: " + ds_args.output_dir + "/dipspades.log\n")

    write_haplocontigs_in_file(ds_args.haplocontigs, ds_args.haplocontigs_fnames)
    write_params(ds_args.output_dir, general_command_line, ds_command_line)

    config_fname = prepare_configs(os.path.join(spades_home, "configs/dipspades"), ds_args, log)
    prepare_config(config_fname, ds_args, log)

    try:
        log.info("===== Assembling started.\n")
        binary_path = check_binary(bin_home, log)
        command = [binary_path, config_fname]
        support.sys_call(command, log)
        log.info("\n===== Assembling finished.\n")
        print_ds_output(ds_args.output_dir, log)
        log.info("\n======= dipSPAdes finished.\n")
        log.info("dipSPAdes log can be found here: " + ds_args.output_dir + "/dipspades.log\n")
        log.info("Thank you for using dipSPAdes!")
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == OSError and exc_value.errno == errno.ENOEXEC: # Exec format error
            support.ds_error("It looks like you are using SPAdes binaries for another platform.\n" + get_spades_binaries_info_message())
        else:
            log.exception(exc_value)
            support.ds_error("exception caught: %s" % exc_type, log)
    except BaseException: # since python 2.5 system-exiting exceptions (e.g. KeyboardInterrupt) are derived from BaseException
        exc_type, exc_value, _ = sys.exc_info()
        log.exception(exc_value)
        support.ds_error("exception caught: %s" % exc_type, log)


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
