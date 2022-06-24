#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from gettext import gettext
from os.path import basename
from os.path import abspath, expanduser

import support
import options_storage
from process_cfg import empty_config

def get_mode():
    mode = None
    script_basename = basename(options_storage.first_command_line[0])
    options = options_storage.first_command_line

    mode_parser = argparse.ArgumentParser(add_help=False)
    mode_parser.add_argument("--isolate", dest="isolate", action="store_true")
    mode_parser.add_argument("--rna", dest="rna", action="store_true")
    mode_parser.add_argument("--plasmid", dest="plasmid", action="store_true")
    mode_parser.add_argument("--meta", dest="meta", action="store_true")
    mode_parser.add_argument("--bio", dest="bio", action="store_true")
    mode_parser.add_argument("--metaviral", dest="metaviral", action="store_true")
    mode_parser.add_argument("--metaplasmid", dest="metaplasmid", action="store_true")
    mode_parser.add_argument("--rnaviral", dest="rnaviral", action="store_true")
    mode_parser.add_argument("--corona", dest="corona", action="store_true")
    nargs, unknown_args = mode_parser.parse_known_args(options)

    if script_basename == "rnaspades.py" or nargs.rna:
        mode = "rna"
    elif script_basename == "rnaviralspades.py" or nargs.rnaviral:
        mode = "rnaviral"
    elif script_basename == "plasmidspades.py" or nargs.plasmid:
        mode = "plasmid"
    elif nargs.bio:
        mode = "bgc"
    elif script_basename == "metaspades.py" or nargs.meta:
        mode = "meta"
    if script_basename == "metaplasmidspades.py" or (nargs.plasmid and nargs.meta) or nargs.metaplasmid:
        mode = "metaplasmid"
    if script_basename == "metaviralspades.py" or nargs.metaviral:
        mode = "metaviral"
    if script_basename == "coronaspades.py" or nargs.corona:
        mode = "corona"
    return mode


def add_mode_to_args(args):
    mode = get_mode()
    if mode == "rna":
        args.rna = True
    elif mode == "plasmid":
        args.plasmid = True
    elif mode == "bgc":
        args.meta = True
        args.bio = True
    elif mode == "meta":
        args.meta = True
    elif mode == "metaplasmid":
        args.meta = True
        args.plasmid = True
        args.metaplasmid = True
    elif mode == "metaviral":
        args.meta = True
        args.plasmid = True
        args.metaviral = True
    elif mode == "rnaviral":
        args.meta = True
        args.rnaviral = True
    elif mode == "corona":
        args.meta = True
        args.rnaviral = True
        args.corona = True



def version():
    mode = get_mode()
    ver = "SPAdes genome assembler v%s" % options_storage.spades_version
    if mode is not None:
        ver += " [%sSPAdes mode]" % mode
    return ver


class SpadesHelpFormatter(argparse.HelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=30, width=100):
        super(SpadesHelpFormatter, self).__init__(prog, indent_increment, max_help_position, width)

    def _split_lines(self, text, width):
        return text.splitlines()

    def _format_usage(self, usage, actions, group, prefix=None):
        if prefix is None:
            prefix = gettext(version() + "\n\nUsage: ")
        return argparse.HelpFormatter._format_usage(self, usage, actions, group, prefix)


def init_dataset_data():
    return dict()


class AddToDatasetAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None,
                 required=False, help=None, metavar=None):
        super(AddToDatasetAction, self).__init__(option_strings, dest, nargs, const, default, type, choices, required,
                                                 help, metavar)

    def __call__(self, parser, namespace, values, option_string=None):
        if option_string == "-s":
            support.old_style_single_reads = True
        if option_string not in options_storage.OLD_STYLE_READS_OPTIONS:
            support.only_old_style_options = False

        # create dataset_data if don't exsist
        if (not "dataset_data" in namespace) or (namespace.dataset_data is None):
            dataset_data = init_dataset_data()
            setattr(namespace, "dataset_data", dataset_data)

        # transfer new format to old
        arg = ""
        if len(values) == 2:
            opt = "--" + option_string.split('-')[2] + values[0]
            if len(option_string.split('-')) > 3:
                if option_string.split('-')[-1] == "or":
                    opt += "-" + values[1]
                else:
                    opt += "-" + option_string.split('-')[-1]
                    arg = values[-1]
            else:
                arg = values[-1]
        else:
            opt = option_string
            if len(values) > 0:
                arg = values[0]

        # add to dataset for old format
        support.add_to_dataset(opt, arg, namespace.dataset_data)


class StoreUniqueAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None,
                 required=False, help=None, metavar=None):
        super(StoreUniqueAction, self).__init__(option_strings=option_strings, dest=dest, nargs=nargs, const=const,
                                                default=default, type=type, choices=choices, required=required,
                                                help=help, metavar=metavar)

    def __call__(self, parser, namespace, values, option_string=None):
        if namespace.__dict__[self.dest] is not None:
            raise argparse.ArgumentError(self, "option was specified at least twice")
        setattr(namespace, self.dest, values)


class ConcatenationAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None,
                 required=False, help=None, metavar=None):
        super(ConcatenationAction, self).__init__(option_strings=option_strings, dest=dest, nargs=nargs, const=const,
                                                  default=default, type=type, choices=choices, required=required,
                                                  help=help, metavar=metavar)

    def __call__(self, parser, namespace, values, option_string=None):
        values = [x for outer in values for x in outer]
        if len(values) == 1 and values[0] == "auto":
            values = values[0]
        elif len(values) > 1 and "auto" in values:
            raise argparse.ArgumentError(self, "cann't set 'auto' and kmers' size at the same time")
        setattr(namespace, self.dest, values)


def kmer(arg):
    if arg == "auto":
        return [arg]
    else:
        k = int(arg)
        if k < options_storage.MIN_K or k > options_storage.MAX_K:
            raise argparse.ArgumentTypeError("wrong k value %d: all k values should be between %d and %d" %
                                             (k, options_storage.MIN_K, options_storage.MAX_K))
        if k % 2 == 0:
            raise argparse.ArgumentTypeError("wrong k value %d: all k values should be odd" % k)
        return [k]


def kmers(arg):
    k_mers = arg
    if k_mers[-1] == ',':
        k_mers = k_mers[:-1]
    k_mers = k_mers.split(",")
    for i in range(len(k_mers)):
        k_mers[i] = kmer(k_mers[i])[0]
    return k_mers


def qvoffset(arg):
    if arg == "auto":
        return arg
    else:
        return int(arg)


def cov_cutoff(arg):
    if arg == "auto" or arg == "off":
        return arg
    elif support.is_float(arg) and float(arg) > 0.0:
        return float(arg)
    else:
        raise argparse.ArgumentTypeError("wrong value %s (should be a positive float number, or 'auto', or 'off')" % arg)


def lcer_cutoff(arg):
    if support.is_float(arg) and float(arg) > 0.0:
        return float(arg)
    else:
        raise argparse.ArgumentTypeError("wrong value %s (should be a positive float number)" % arg)


def restart_from(arg):
    if arg not in options_storage.SHORT_STAGES_NAME and arg != options_storage.LAST_STAGE and not arg.startswith("k"):
        raise argparse.ArgumentTypeError("wrong value %s (should be 'ec', 'as', 'k<int>', 'mc', or '%s')" % (arg, options_storage.LAST_STAGE))
    return arg


def stop_after(arg):
    if arg not in options_storage.SHORT_STAGES_NAME and not arg.startswith("k"):
        raise argparse.ArgumentTypeError("wrong value %s (should be 'ec', 'as', 'k<int>', or 'mc')" % arg)
    return arg


def read_cov_threshold(arg):
    if support.is_int(arg) and int(arg) >= 0:
        return int(arg)
    else:
        raise argparse.ArgumentTypeError("wrong value %s (should be a non-negative integer number)" % arg)


def add_deprecated_input_data_args(pgroup_input_data):
    for num in range(1, 10):
        for sufix in ["-12", "-1", "-2", "-s"]:
            pgroup_input_data.add_argument("--pe%d%s" % (num, sufix),
                                           metavar="<filename>",
                                           nargs=1,
                                           help=argparse.SUPPRESS,
                                           action=AddToDatasetAction)
            pgroup_input_data.add_argument("--mp%d%s" % (num, sufix),
                                           metavar="<filename>",
                                           nargs=1,
                                           help=argparse.SUPPRESS,
                                           action=AddToDatasetAction)
            pgroup_input_data.add_argument("--hqmp%d%s" % (num, sufix),
                                           metavar="<filename>",
                                           nargs=1,
                                           help=argparse.SUPPRESS,
                                           action=AddToDatasetAction)

        for orientation in ["-fr", "-rf", "-ff"]:
            pgroup_input_data.add_argument("--pe%d%s" % (num, orientation),
                                           nargs=0,
                                           help=argparse.SUPPRESS,
                                           action=AddToDatasetAction)
            pgroup_input_data.add_argument("--mp%d%s" % (num, orientation),
                                           nargs=0,
                                           help=argparse.SUPPRESS,
                                           action=AddToDatasetAction)
            pgroup_input_data.add_argument("--hqmp%d%s" % (num, orientation),
                                           nargs=0,
                                           help=argparse.SUPPRESS,
                                           action=AddToDatasetAction)

        pgroup_input_data.add_argument("--s%d" % num,
                                       metavar="<filename>",
                                       nargs=1,
                                       help=argparse.SUPPRESS,
                                       action=AddToDatasetAction)
        pgroup_input_data.add_argument("--pe%d-m" % num,
                                       metavar="<filename>",
                                       nargs=1,
                                       help=argparse.SUPPRESS,
                                       action=AddToDatasetAction)


def add_basic_args(pgroup_basic):
    mode = get_mode()
    pgroup_basic.add_argument("-o",
                              metavar="<output_dir>",
                              help="directory to store all the resulting files (required)",
                              type=str,
                              default=None,
                              dest="output_dir",
                              action=StoreUniqueAction)

    help_hidden = (mode is not None)
    pgroup_basic.add_argument("--isolate",
                              dest="isolate",
                              help="this flag is highly recommended for high-coverage isolate and multi-cell data"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--sc",
                              dest="single_cell",
                              help="this flag is required for MDA (single-cell) data"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--meta",
                              dest="meta",
                              help="this flag is required for metagenomic data"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--bio",
                              dest="bio",
                              help="this flag is required for biosyntheticSPAdes mode"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--corona",
                              dest="corona",
                              help="this flag is required for coronaSPAdes mode"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--rna",
                              dest="rna",
                              help="this flag is required for RNA-Seq data"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--plasmid",
                              dest="plasmid",
                              help="runs plasmidSPAdes pipeline for plasmid detection"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--metaviral",
                              dest="metaviral",
                              help="runs metaviralSPAdes pipeline for virus detection"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--metaplasmid",
                              dest="metaplasmid",
                              help="runs metaplasmidSPAdes pipeline for plasmid detection in metagenomic datasets (equivalent for --meta --plasmid)"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")                              
    pgroup_basic.add_argument("--rnaviral",
                              dest="rnaviral",
                              help="this flag enables virus assembly module from RNA-Seq data"
                              if not help_hidden else argparse.SUPPRESS,
                              action="store_true")
    pgroup_basic.add_argument("--iontorrent",
                              dest="iontorrent",
                              help="this flag is required for IonTorrent data",
                              action="store_true")
    pgroup_basic.add_argument("--test",
                              dest="test_mode",
                              help="runs SPAdes on toy dataset",
                              action="store_true")
    pgroup_basic.add_argument("-h", "--help",
                              help="prints this usage message",
                              action="help")
    pgroup_basic.add_argument("-v", "--version",
                              help="prints version",
                              action="version",
                              version=version())


def add_library_args(libid, name, suffixes, pgroup_input_data, help_hidden=False):
    if "12" in suffixes:
        pgroup_input_data.add_argument("--%s-12" % libid,
                                       metavar=("<#>", "<filename>"),
                                       nargs=2,
                                       help="file with interlaced reads for %s library number <#>.\n"
                                            "Older deprecated syntax is -%s<#>-12 <filename>" % (name, libid)
                                       if not help_hidden else argparse.SUPPRESS,
                                       action=AddToDatasetAction)
    if "1" in suffixes:
        pgroup_input_data.add_argument("--%s-1" % libid, metavar=("<#>", "<filename>"),
                                       nargs=2,
                                       help="file with forward reads for %s library number <#>.\n"
                                            "Older deprecated syntax is -%s<#>-1 <filename>" % (name, libid)
                                       if not help_hidden else argparse.SUPPRESS,
                                       action=AddToDatasetAction)

    if "2" in suffixes:
        pgroup_input_data.add_argument("--%s-2" % libid,
                                       metavar=("<#>", "<filename>"),
                                       nargs=2,
                                       help="file with reverse reads for %s library number <#>.\n"
                                            "Older deprecated syntax is -%s<#>-2 <filename>" % (name, libid)
                                       if not help_hidden else argparse.SUPPRESS,
                                       action=AddToDatasetAction)

    if "s" in suffixes:
        pgroup_input_data.add_argument("--%s-s" % libid,
                                       metavar=("<#>", "<filename>"),
                                       nargs=2,
                                       help="file with unpaired reads for %s library number <#>.\n"
                                            "Older deprecated syntax is -%s<#>-s <filename>" % (name, libid)
                                       if not help_hidden else argparse.SUPPRESS,
                                       action=AddToDatasetAction)

    if "m" in suffixes:
        pgroup_input_data.add_argument("--%s-m" % libid,
                                       metavar=("<#>", "<filename>"),
                                       nargs=2,
                                       help="file with merged reads for %s library number <#>.\n"
                                            "Older deprecated syntax is -%s<#>-m <filename>" % (name, libid)
                                       if not help_hidden else argparse.SUPPRESS,
                                       action=AddToDatasetAction)

    if "or" in suffixes:
        pgroup_input_data.add_argument("--%s-or" % libid,
                                       metavar=("<#>", "<or>"),
                                       nargs=2,
                                       help="orientation of reads for %s library number <#> \n(<or> = fr, rf, ff).\n"
                                            "Older deprecated syntax is -%s<#>-<or>" % (name, libid)
                                       if not help_hidden else argparse.SUPPRESS,
                                       action=AddToDatasetAction)


def add_input_data_args(pgroup_input_data):
    mode = get_mode()

    pgroup_input_data.add_argument("--12",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with interlaced forward and reverse paired-end reads",
                                   action=AddToDatasetAction)
    pgroup_input_data.add_argument("-1",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with forward paired-end reads",
                                   action=AddToDatasetAction)
    pgroup_input_data.add_argument("-2",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with reverse paired-end reads",
                                   action=AddToDatasetAction)
    pgroup_input_data.add_argument("-s",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with unpaired reads",
                                   action=AddToDatasetAction)
    pgroup_input_data.add_argument("--merged",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with merged forward and reverse paired-end reads",
                                   action=AddToDatasetAction)

    add_deprecated_input_data_args(pgroup_input_data)
    help_hidden = (mode in ["rna", "meta"])
    add_library_args("pe", "paired-end", ["12", "1", "2", "s", "m", "or"], pgroup_input_data)
    pgroup_input_data.add_argument("--s",
                                   metavar=("<#>", "<filename>"),
                                   nargs=2,
                                   help="file with unpaired reads for single reads library number <#>.\n"
                                        "Older deprecated syntax is --s<#> <filename>",
                                   action=AddToDatasetAction)
    add_library_args("mp", "mate-pair", ["12", "1", "2", "s", "or"], pgroup_input_data, help_hidden)
    add_library_args("hqmp", "high-quality mate-pair", ["12", "1", "2", "s", "or"], pgroup_input_data, help_hidden)

    pgroup_input_data.add_argument("--sanger",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with Sanger reads"
                                   if not help_hidden else argparse.SUPPRESS,
                                   action=AddToDatasetAction)

    help_hidden = (mode == "rna")
    pgroup_input_data.add_argument("--pacbio",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with PacBio reads",
                                   action=AddToDatasetAction)
    pgroup_input_data.add_argument("--nanopore",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with Nanopore reads",
                                   action=AddToDatasetAction)

    help_hidden = (mode == "meta")
    pgroup_input_data.add_argument("--trusted-contigs",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with trusted contigs"
                                   if not help_hidden else argparse.SUPPRESS,
                                   action=AddToDatasetAction)

    pgroup_input_data.add_argument("--untrusted-contigs",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with untrusted contigs"
                                   if not help_hidden else argparse.SUPPRESS,
                                   action=AddToDatasetAction)

    help_hidden = (mode != "rna")
    pgroup_input_data.add_argument("--fl-rna",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with PacBio/Nanopore/contigs that capture full-length transcripts"
                                   if not help_hidden else argparse.SUPPRESS,
                                   action=AddToDatasetAction)
    pgroup_input_data.add_argument("--ss",
                                   metavar="<type>",
                                   dest="strand_specificity",
                                   choices=["fr", "rf"],
                                   help="strand specific data, <type> = fr (normal) and rf (antisense).\n"
                                        "Older deprecated syntax is --ss-<type>"
                                   if not help_hidden else argparse.SUPPRESS,
                                   action="store")

    pgroup_input_data.add_argument("--ss-fr",
                                   metavar="<type>",
                                   dest="strand_specificity",
                                   const="fr",
                                   help=argparse.SUPPRESS,
                                   action="store_const")

    pgroup_input_data.add_argument("--ss-rf",
                                   dest="strand_specificity",
                                   const="rf",
                                   help=argparse.SUPPRESS,
                                   action="store_const")

    help_hidden = (mode != "plasmid" and mode != "bgc" and mode != "metaplasmid" and mode != "metaviral")
    pgroup_input_data.add_argument("--assembly-graph",
                                   metavar="<filename>",
                                   nargs=1,
                                   help="file with assembly graph"
                                   if not help_hidden else argparse.SUPPRESS,
                                   action=AddToDatasetAction)
    
def add_pipeline_args(pgroup_pipeline):
    mode = get_mode()
    help_hidden = (mode == "rna" or mode == "rnaviral")
    pgroup_pipeline.add_argument("--only-error-correction",
                                 dest="only_error_correction",
                                 default=None,
                                 help="runs only read error correction (without assembling)"
                                 if not help_hidden else argparse.SUPPRESS,
                                 action="store_true")
    pgroup_pipeline.add_argument("--only-assembler",
                                 dest="only_assembler",
                                 default=None,
                                 help="runs only assembling (without read error correction)"
                                 if not help_hidden else argparse.SUPPRESS,
                                 action="store_true")

    help_hidden = (mode in ["rna", "meta", "rnaviral"])
    careful_group = pgroup_pipeline.add_mutually_exclusive_group()
    careful_group.add_argument("--careful",
                               dest="careful",
                               default=None,
                               help="tries to reduce number of mismatches and short indels"
                               if not help_hidden else argparse.SUPPRESS,
                               action="store_true")
    careful_group.add_argument("--careful:false",
                               dest="careful",
                               default=None,
                               help=argparse.SUPPRESS,
                               action="store_false")
    pgroup_pipeline.add_argument("--checkpoints",
                                 metavar="<last or all>",
                                 dest="checkpoints",
                                 help="save intermediate check-points ('last', 'all')",
                                 action="store")
    pgroup_pipeline.add_argument("--continue",
                                 dest="continue_mode",
                                 help="continue run from the last available check-point (only -o should be specified)",
                                 action="store_true")

    restart_from_help = "restart run with updated options and from the specified check-point\n" \
                        "('ec', 'as', 'k<int>', 'mc', '%s')" % options_storage.LAST_STAGE
    if mode in ["rna", "rnaviral"]:
        restart_from_help = "restart run with updated options and from the specified check-point\n" \
                            "('as', 'k<int>', '%s')" % options_storage.LAST_STAGE
    pgroup_pipeline.add_argument("--restart-from",
                                 metavar="<cp>",
                                 dest="restart_from",
                                 default=None,
                                 type=restart_from,
                                 help=restart_from_help,
                                 action="store")

    disable_gzip_output_group = pgroup_pipeline.add_mutually_exclusive_group()
    disable_gzip_output_group.add_argument("--disable-gzip-output",
                                           dest="disable_gzip_output",
                                           default=None,
                                           help="forces error correction not to compress the corrected reads",
                                           action="store_true")
    disable_gzip_output_group.add_argument("--disable-gzip-output:false",
                                           dest="disable_gzip_output",
                                           default=None,
                                           help=argparse.SUPPRESS,
                                           action="store_false")

    disable_rr = pgroup_pipeline.add_mutually_exclusive_group()
    disable_rr.add_argument("--disable-rr",
                            dest="disable_rr",
                            default=None,
                            help="disables repeat resolution stage of assembling",
                            action="store_true")
    disable_rr.add_argument("--disable-rr:false",
                            dest="disable_rr",
                            default=None,
                            help=argparse.SUPPRESS,
                            action="store_false")


def add_advanced_args(pgroup_advanced):
    mode = get_mode()
    pgroup_advanced.add_argument("--dataset",
                                 metavar="<filename>",
                                 type=support.check_file_existence,
                                 dest="dataset_yaml_filename",
                                 help="file with dataset description in YAML format",
                                 action="store")

    pgroup_advanced.add_argument("-t", "--threads",
                                 metavar="<int>",
                                 dest="threads",
                                 type=int,
                                 help="number of threads. [default: %s]\n" % options_storage.THREADS,
                                 action="store")

    pgroup_advanced.add_argument("-m", "--memory",
                                 metavar="<int>",
                                 type=int,
                                 dest="memory",
                                 help="RAM limit for SPAdes in Gb (terminates if exceeded). [default: %s]\n" % options_storage.MEMORY,
                                 action="store")
    pgroup_advanced.add_argument("--tmp-dir",
                                 metavar="<dirname>",
                                 help="directory for temporary files. [default: <output_dir>/tmp]",
                                 dest="tmp_dir",
                                 action="store")

    pgroup_advanced.add_argument("-k",
                                 metavar="<int>",
                                 dest="k_mers",
                                 nargs='+',
                                 type=kmers,
                                 help="list of k-mer sizes (must be odd and less than %d)\n"
                                      "[default: 'auto']" % (options_storage.MAX_K + 1),
                                 action=ConcatenationAction)

    help_hidden = (mode in ["rna", "meta", "rnaviral"])
    pgroup_advanced.add_argument("--cov-cutoff",
                                 metavar="<float>",
                                 type=cov_cutoff,
                                 default=None,
                                 dest="cov_cutoff",
                                 help="coverage cutoff value (a positive float number, "
                                      "or 'auto', or 'off')\n[default: 'off']"
                                 if not help_hidden else argparse.SUPPRESS,
                                 action="store")

    pgroup_advanced.add_argument("--phred-offset",
                                 metavar="<33 or 64>",
                                 dest="qvoffset",
                                 type=qvoffset,
                                 help="PHRED quality offset in the input reads (33 or 64),\n"
                                      "[default: auto-detect]",
                                 action="store")

    pgroup_advanced.add_argument("--custom-hmms",
                                 metavar="<dirname>",
                                 dest="custom_hmms",
                                 help="directory with custom hmms that replace default ones,\n"
                                      "[default: None]",
                                 action="store")


def add_hidden_args(pgroup_hidden):
    show_help_hidden = ("--help-hidden" in sys.argv)

    debug_group = pgroup_hidden.add_mutually_exclusive_group()
    debug_group.add_argument("--debug",
                             dest="developer_mode",
                             default=None,
                             help="runs SPAdes in debug mode"
                             if show_help_hidden else argparse.SUPPRESS,
                             action="store_true")
    debug_group.add_argument("--debug:false",
                             dest="developer_mode",
                             default=None,
                             help=argparse.SUPPRESS,
                             action="store_false")
    debug_group.add_argument("--trace-time",
                             dest="time_tracer",
                             default=None,
                             help="enable time tracker"
                             if show_help_hidden else argparse.SUPPRESS,
                             action="store_true")

    pgroup_hidden.add_argument("--stop-after",
                               metavar="<cp>",
                               dest="stop_after",
                               type=stop_after,
                               help="runs SPAdes until the specified check-point ('ec', 'as', 'k<int>', 'mc') inclusive"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store")
    pgroup_hidden.add_argument("--truseq",
                               dest="truseq_mode",
                               default=None,
                               help="runs SPAdes in TruSeq mode"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store_true")

    mismatch_correction_group = pgroup_hidden.add_mutually_exclusive_group()
    mismatch_correction_group.add_argument("--mismatch-correction",
                                           dest="mismatch_corrector",
                                           default=None,
                                           help="runs post processing correction of mismatches and short indels"
                                           if show_help_hidden else argparse.SUPPRESS,
                                           action="store_true")
    mismatch_correction_group.add_argument("--mismatch-correction:false",
                                           dest="mismatch_corrector",
                                           default=None,
                                           help=argparse.SUPPRESS,
                                           action="store_false")

    pgroup_hidden.add_argument("--reference",
                               metavar="<filename>",
                               dest="reference",
                               type=support.check_file_existence,
                               help="file with reference for deep analysis (only in debug mode)"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store")
    pgroup_hidden.add_argument("--series-analysis",
                               metavar="<filename>",
                               dest="series_analysis",
                               type=support.check_file_existence,
                               help="config for metagenomics-series-augmented reassembly"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store")
    pgroup_hidden.add_argument("--configs-dir",
                               metavar="<config_dir>",
                               dest="configs_dir",
                               type=support.check_dir_existence,
                               help="directory with configs"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store")
    pgroup_hidden.add_argument("--read-buffer-size",
                               metavar="<int>",
                               dest="read_buffer_size",
                               type=int,
                               help="sets size of read buffer for graph construction"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store")
    pgroup_hidden.add_argument("--large-genome",
                               dest="large_genome",
                               default=False,
                               help="Enables optimizations for large genomes"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store_true")
    pgroup_hidden.add_argument("--save-gp",
                               dest="save_gp",
                               default=None,
                               help="Enables saving graph pack before repeat resolution (even without --debug)"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store_true")
    pgroup_hidden.add_argument("--hidden-cov-cutoff",
                               metavar="<float>",
                               type=lcer_cutoff,
                               dest="lcer_cutoff",
                               help="coverage cutoff value deeply integrated in simplification" \
                                    " (a positive float number). Base coverage! Will be adjusted depending on K and RL!"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store")
    pgroup_hidden.add_argument("--read-cov-threshold",
                               metavar="<int>",
                               dest="read_cov_threshold",
                               type=read_cov_threshold,
                               help="read median coverage threshold (non-negative integer)"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store")
    pgroup_hidden.add_argument("--only-generate-config",
                               dest="only_generate_config",
                               help="generate configs and print script to run_spades.sh"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="store_true")
    pgroup_hidden.add_argument("--no-clear-after",
                               dest="no_clear_after",
                               help="don't delete tmp files after SPAdes pipeline finished"
                               if show_help_hidden else argparse.SUPPRESS,
                               action = "store_true")
    pgroup_hidden.add_argument("--help-hidden",
                               help="prints this usage message with all hidden options"
                               if show_help_hidden else argparse.SUPPRESS,
                               action="help")


def create_parser():
    parser = argparse.ArgumentParser(prog="spades.py", formatter_class=SpadesHelpFormatter,
                                     usage="%(prog)s [options] -o <output_dir>", add_help=False)

    #pgroup for parser group
    pgroup_basic = parser.add_argument_group("Basic options")
    pgroup_input_data = parser.add_argument_group("Input data")
    pgroup_pipeline = parser.add_argument_group("Pipeline options")
    pgroup_advanced = parser.add_argument_group("Advanced options")
    pgroup_hidden = parser.add_argument_group("Hidden options")

    add_basic_args(pgroup_basic)
    add_input_data_args(pgroup_input_data)
    add_pipeline_args(pgroup_pipeline)
    add_advanced_args(pgroup_advanced)
    add_hidden_args(pgroup_hidden)

    return parser


def check_options_for_restart_from(log):
    if ("dataset_data" in options_storage.args) and (options_storage.args.dataset_data is not None):
        support.error("you cannot specify input data (-1, -2, -12, --pe-1, --pe-2 ...) with --restart-from option!", log)
    if options_storage.args.dataset_yaml_filename:
        support.error("you cannot specify --dataset with --restart-from option!", log)
    if options_storage.args.single_cell:
        support.error("you cannot specify --sc with --restart-from option!", log)
    if options_storage.args.meta:
        support.error("you cannot specify --meta with --restart-from option!", log)
    if options_storage.args.plasmid:
        support.error("you cannot specify --plasmid with --restart-from option!", log)
    if options_storage.args.rna:
        support.error("you cannot specify --rna with --restart-from option!", log)
    if options_storage.args.isolate:
        support.error("you cannot specify --isolate with --restart-from option!", log)
    if options_storage.args.iontorrent:
        support.error("you cannot specify --iontorrent with --restart-from option!", log)
    if options_storage.args.only_assembler:
        support.error("you cannot specify --only-assembler with --restart-from option!", log)
    if options_storage.args.only_error_correction:
        support.error("you cannot specify --only-error-correction with --restart-from option!", log)
    if options_storage.args.strand_specificity is not None:
        support.error("you cannot specify strand specificity (--ss-rf or --ss-fr) with --restart-from option!", log)

def add_to_option(args, log, skip_output_dir):
    if args.restart_from:
        check_options_for_restart_from(log)

    add_mode_to_args(options_storage.args)

    if args.test_mode:
        if not skip_output_dir:
            if "output_dir" in options_storage.args and options_storage.args.output_dir is not None:
                support.error("you cannot specify -o and --test simultaneously")
            options_storage.args.output_dir = os.path.abspath("spades_test")
    else:
        if "output_dir" not in options_storage.args or options_storage.args.output_dir is None:
            support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).")

    if not skip_output_dir:
        output_dir = abspath(expanduser(args.output_dir))
        options_storage.dict_of_rel2abs[args.output_dir] = output_dir
        support.check_path_is_ascii(output_dir, "output directory")
        args.output_dir = output_dir

    if args.tmp_dir is not None:
        tmp_dir = abspath(expanduser(args.tmp_dir))
        options_storage.dict_of_rel2abs[args.tmp_dir] = tmp_dir
        support.check_path_is_ascii(tmp_dir, "directory for temporary files")
        args.tmp_dir = tmp_dir

    if args.custom_hmms is not None:
        custom_hmms = abspath(expanduser(args.custom_hmms))
        support.check_path_is_ascii(custom_hmms, "directory with custom hmms")
        args.custom_hmms = custom_hmms

    if "reference" in args and args.reference is not None:
        args.developer_mode = True

    if args.developer_mode:
        args.no_clear_after = True

    if args.only_assembler and args.only_error_correction:
        support.error("you cannot specify --only-error-correction and --only-assembler simultaneously")

    if (args.rna or args.rnaviral) and args.only_error_correction:
        support.error("you cannot specify --only-error-correction in RNA-seq mode!", log)

    if args.isolate and args.only_error_correction:
        support.error("you cannot specify --only-error-correction in isolate mode!", log)

    if args.careful == False and args.mismatch_corrector == True:
        support.error("you cannot specify --mismatch-correction and --careful:false simultaneously")

    if args.careful == True and args.mismatch_corrector == False:
        support.error("you cannot specify --mismatch-correction:false and --careful simultaneously")

    if args.rna and (args.careful or args.mismatch_corrector):
        support.error("you cannot specify --mismatch-correction or --careful in RNA-seq mode!", log)

    if args.isolate and (args.careful or args.mismatch_corrector):
        support.error("you cannot specify --mismatch-correction or --careful in isolate mode!", log)

    if args.only_assembler and args.isolate:
        support.warning("Isolate mode already implies --only-assembler, so this option has no effect.")

    if args.only_assembler and args.rna:
        support.warning("RNA mode already implies --only-assembler, so this option has no effect.")

    if args.restart_from is not None:
        args.continue_mode = True
    if args.careful is not None:
        args.mismatch_corrector = args.careful
    if args.truseq_mode:
        enable_truseq_mode()
    if (args.isolate or args.rna or args.rnaviral) and not args.iontorrent:
        args.only_assembler = True


def add_to_cfg(cfg, log, bin_home, spades_home, args):
    ### FILLING cfg
    cfg["common"] = empty_config()
    cfg["dataset"] = empty_config()
    if not args.only_assembler:
        cfg["error_correction"] = empty_config()
    if not args.only_error_correction:
        cfg["assembly"] = empty_config()

    # common
    cfg["common"].__dict__["checkpoints"] = args.checkpoints
    cfg["common"].__dict__["output_dir"] = args.output_dir
    cfg["common"].__dict__["tmp_dir"] = args.tmp_dir
    cfg["common"].__dict__["max_threads"] = args.threads
    cfg["common"].__dict__["max_memory"] = args.memory
    cfg["common"].__dict__["developer_mode"] = args.developer_mode
    cfg["common"].__dict__["time_tracer"] = args.time_tracer
    if args.series_analysis:
        cfg["common"].__dict__["series_analysis"] = args.series_analysis

    hmms_path = None
    if args.custom_hmms is not None:
        hmms_path = args.custom_hmms
    elif args.bio:
        hmms_path = os.path.join(spades_home, options_storage.biosyntheticspades_hmms)
    elif args.corona:
        hmms_path = os.path.join(spades_home, options_storage.coronaspades_hmms)

    if hmms_path is not None:
        hmms = ""
        is_hmmfile= lambda hmmfile: os.path.isfile(hmmfile) \
            and (hmmfile.endswith("hmm") or hmmfile.endswith("hmm.gz") or \
                 hmmfile.endswith("aa") or hmmfile.endswith("aa.gz") or \
                 hmmfile.endswith("fa") or hmmfile.endswith("fa.gz") or \
                 hmmfile.endswith("fna") or hmmfile.endswith("fna.gz"))
        if os.path.isdir(hmms_path):
            hmms = ",".join([os.path.join(hmms_path, hmmfile)
                             for hmmfile in os.listdir(hmms_path)
                             if is_hmmfile(os.path.join(hmms_path, hmmfile))])
        elif is_hmmfile(hmms_path):
            hmms = hmms_path

        if hmms == "":
            support.error("Custom HMM folder does not contain any HMMs. They should have .hmm or .hmm.gz extension.", log)
        cfg["common"].__dict__["set_of_hmms"] = hmms

    # dataset section
    cfg["dataset"].__dict__["yaml_filename"] = args.dataset_yaml_filename
    if args.developer_mode and args.reference:
        cfg["dataset"].__dict__["reference"] = args.reference

    # error correction
    if not args.only_assembler:
        cfg["error_correction"].__dict__["output_dir"] = os.path.join(cfg["common"].output_dir, "corrected")
        cfg["error_correction"].__dict__["gzip_output"] = not args.disable_gzip_output
        if args.qvoffset:
            cfg["error_correction"].__dict__["qvoffset"] = args.qvoffset
        cfg["error_correction"].__dict__["iontorrent"] = args.iontorrent
        cfg["error_correction"].__dict__["max_iterations"] = options_storage.ITERATIONS
        if args.meta or args.large_genome:
            cfg["error_correction"].__dict__["count_filter_singletons"] = 1
        if args.read_buffer_size:
            cfg["error_correction"].__dict__["read_buffer_size"] = args.read_buffer_size

    # assembly
    if not args.only_error_correction:
        if args.k_mers == "auto" and args.restart_from is None:
            args.k_mers = None
        if args.k_mers:
            cfg["assembly"].__dict__["iterative_K"] = args.k_mers
        elif (args.rna or args.rnaviral):
            cfg["assembly"].__dict__["iterative_K"] = "auto"
        else:
            cfg["assembly"].__dict__["iterative_K"] = options_storage.K_MERS_SHORT
        cfg["assembly"].__dict__["disable_rr"] = args.disable_rr
        cfg["assembly"].__dict__["cov_cutoff"] = args.cov_cutoff
        cfg["assembly"].__dict__["lcer_cutoff"] = args.lcer_cutoff
        cfg["assembly"].__dict__["save_gp"] = args.save_gp
        if args.read_buffer_size:
            cfg["assembly"].__dict__["read_buffer_size"] = args.read_buffer_size
        cfg["assembly"].__dict__["correct_scaffolds"] = options_storage.correct_scaffolds

    # corrector can work only if contigs exist (not only error correction)
    if (not args.only_error_correction) and args.mismatch_corrector:
        cfg["mismatch_corrector"] = empty_config()
        cfg["mismatch_corrector"].__dict__["skip-masked"] = None
        cfg["mismatch_corrector"].__dict__["bwa"] = os.path.join(bin_home, "spades-bwa")
        cfg["mismatch_corrector"].__dict__["threads"] = args.threads
        cfg["mismatch_corrector"].__dict__["output-dir"] = args.output_dir
    cfg["run_truseq_postprocessing"] = options_storage.run_truseq_postprocessing


def postprocessing(args, cfg, dataset_data, log, spades_home, load_processed_dataset, restart_from, options=None):
    if sys.version.startswith("2."):
        import pyyaml2 as pyyaml
    elif sys.version.startswith("3."):
        import pyyaml3 as pyyaml

    if args.test_mode:
        if args.plasmid:
            support.add_to_dataset("-1", os.path.join(spades_home, "test_dataset_plasmid/pl1.fq.gz"), dataset_data)
            support.add_to_dataset("-2", os.path.join(spades_home, "test_dataset_plasmid/pl2.fq.gz"), dataset_data)
        else:
            support.add_to_dataset("-1", os.path.join(spades_home, "test_dataset/ecoli_1K_1.fq.gz"), dataset_data)
            support.add_to_dataset("-2", os.path.join(spades_home, "test_dataset/ecoli_1K_2.fq.gz"), dataset_data)

    if args.bio or args.rnaviral:
        args.meta = True
    if not args.output_dir:
        support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).", log)
    if not os.path.isdir(args.output_dir):
        if args.continue_mode:
            support.error("the output_dir should exist for --continue and for --restart-from!", log)
        os.makedirs(args.output_dir)
    if args.restart_from or restart_from:
        if args.continue_mode:  # saving parameters specified with --restart-from
            if not support.dataset_is_empty(dataset_data):
                support.error("you cannot specify reads with --restart-from option!", log)
            save_restart_options()
        else:  # overriding previous run parameters
            load_restart_options()
    elif args.continue_mode:  # it is just --continue, NOT --restart-from
        continue_parser = argparse.ArgumentParser(add_help=False)
        continue_parser.add_argument("--continue", dest="continue_mode", action="store_true")
        continue_parser.add_argument("-o", type=str, dest="output_dir", action=StoreUniqueAction)
        nargs, unknown_args = continue_parser.parse_known_args(options)
        if unknown_args:
            support.error("you cannot specify any option except -o with --continue option! "
                          "Please use '--restart-from last' if you need to change some "
                          "of the options from the initial run and continue from the last available checkpoint.", log)
    if args.meta:
        if args.careful or args.mismatch_corrector or (args.cov_cutoff != "off" and args.cov_cutoff is not None):
            support.error("you cannot specify --careful, --mismatch-correction or --cov-cutoff in metagenomic mode!",
                          log)
    if args.rna or args.rnaviral:
        if args.careful:
            support.error("you cannot specify --careful in RNA-Seq mode!", log)

    modes_count =  [args.large_genome, args.rna, args.plasmid, args.meta, args.single_cell, args.isolate, args.rnaviral,
                    args.corona, args.metaviral, args.metaplasmid, args.bio].count(True)
    is_metaplasmid = modes_count == 3 and [args.meta, args.plasmid, args.metaplasmid].count(True) == 3
    is_bgc = modes_count == 2 and [args.meta, args.bio].count(True) == 2
    is_metaviral = modes_count == 3 and [args.meta, args.plasmid, args.metaviral].count(True) == 3
    is_rnaviral = modes_count == 2 and [args.meta, args.rnaviral].count(True) == 2
    is_corona = modes_count == 3 and [args.meta, args.rnaviral, args.corona].count(True) == 3
    if not (modes_count <= 1 or is_metaplasmid or is_bgc or is_metaviral or is_rnaviral or is_corona):
        # correct cases:
        # - either there is 1 or 0 modes specified
        # - or there is 1 of 5 allowed combinations
        # everything else is forbidden
        support.error("Specified mode combination is not supported! Check out user manual for available modes.", log)
    elif modes_count == 0:
        support.warning("No assembly mode was specified! If you intend to assemble high-coverage multi-cell/isolate data, use '--isolate' option.")

    if args.continue_mode:
        return None

    existing_dataset_data = None
    processed_dataset_fpath = os.path.join(args.output_dir, "input_dataset.yaml")
    if load_processed_dataset:
        if os.path.isfile(processed_dataset_fpath):
            try:
                existing_dataset_data = pyyaml.load(open(processed_dataset_fpath))
            except pyyaml.YAMLError:
                existing_dataset_data = None

    options_storage.original_dataset_data = dataset_data
    if args.dataset_yaml_filename:
        try:
            options_storage.original_dataset_data = pyyaml.load(open(args.dataset_yaml_filename))
        except pyyaml.YAMLError:
            _, exc, _ = sys.exc_info()
            support.error(
                    "exception caught while parsing YAML file (%s):\n" % args.dataset_yaml_filename + str(exc))
        options_storage.original_dataset_data = support.relative2abs_paths(options_storage.original_dataset_data,
                                                      os.path.dirname(args.dataset_yaml_filename))
    else:
        options_storage.original_dataset_data = support.correct_dataset(options_storage.original_dataset_data)
        options_storage.original_dataset_data = support.relative2abs_paths(options_storage.original_dataset_data, os.getcwd())

    if existing_dataset_data is not None:
        dataset_data = existing_dataset_data
    else:
        dataset_data = options_storage.original_dataset_data

    args.dataset_yaml_filename = processed_dataset_fpath

    support.check_dataset_reads(dataset_data, (args.only_assembler or args.rna), args.iontorrent, log)
    if not support.get_lib_ids_by_type(dataset_data, options_storage.READS_TYPES_USED_IN_CONSTRUCTION):
        support.error("you should specify at least one unpaired, paired-end, or high-quality mate-pairs library!")
    if args.rna:
        if len(dataset_data) != len(
                support.get_lib_ids_by_type(dataset_data, options_storage.READS_TYPES_USED_IN_RNA_SEQ)):
            support.error("you cannot specify any data types except " +
                          ", ".join(options_storage.READS_TYPES_USED_IN_RNA_SEQ) + " in RNA-Seq mode!")
            # if len(support.get_lib_ids_by_type(dataset_data, 'paired-end')) > 1:
            #    support.error('you cannot specify more than one paired-end library in RNA-Seq mode!')
    if args.meta and not args.only_error_correction and not args.rnaviral:
        paired_end_libs = max(1, len(support.get_lib_ids_by_type(dataset_data, "paired-end")))
        graph_libs = max(1, len(support.get_lib_ids_by_type(dataset_data, "assembly-graph")))
        long_read_libs = max(1, len(support.get_lib_ids_by_type(dataset_data, ["pacbio", "nanopore"])))
        
        if len(dataset_data) > paired_end_libs + graph_libs + long_read_libs:
            support.error("you cannot specify any data types except a single paired-end library "
                          "(optionally accompanied by a single library of "
                          "PacBio reads or Nanopore reads) in metaSPAdes mode!")

    if existing_dataset_data is None:
        with open(args.dataset_yaml_filename, 'w') as f:
            pyyaml.dump(dataset_data, f,
                        default_flow_style=False, default_style='"', width=float("inf"))

    set_default_values()
    return dataset_data


def parse_args(log, bin_home, spades_home, secondary_filling, restart_from=False, options=None):
    cfg = dict()
    parser = create_parser()

    if secondary_filling:
        old_output_dir = options_storage.args.output_dir
        old_stop_after = options_storage.args.stop_after

    skip_output_dir = secondary_filling
    load_processed_dataset = secondary_filling

    options_storage.args, argv = parser.parse_known_args(options)

    if options_storage.args.restart_from is not None and not secondary_filling:
        for arg in options_storage.args.__dict__:
            parser.set_defaults(**{arg: None})
        options_storage.args, argv = parser.parse_known_args(options)

    if argv:
        msg = "Please specify option (e.g. -1, -2, -s, etc)) for the following paths: %s"
        parser.error(msg % ", ".join(argv))

    if secondary_filling:
        options_storage.args.output_dir = old_output_dir
        options_storage.args.stop_after = old_stop_after

    add_to_option(options_storage.args, log, skip_output_dir)

    if "dataset_data" in options_storage.args:
        dataset_data = options_storage.args.dataset_data
    else:
        dataset_data = init_dataset_data()
    dataset_data = postprocessing(options_storage.args, cfg, dataset_data, log, spades_home,
                                  load_processed_dataset, restart_from, options)

    if options_storage.args.continue_mode:
        return options_storage.args, None, None

    add_to_cfg(cfg, log, bin_home, spades_home, options_storage.args)
    return options_storage.args, cfg, dataset_data


def usage(spades_version, show_hidden=False, mode=None):
    parser = create_parser()
    parser.print_help()


def set_default_values():
    if options_storage.args.threads is None:
        options_storage.args.threads = options_storage.THREADS
    if options_storage.args.memory is None:
        if support.get_available_memory():
            options_storage.args.memory = int(min(options_storage.MEMORY, support.get_available_memory()))
        else:
            options_storage.args.memory = options_storage.MEMORY
    if options_storage.args.disable_gzip_output is None:
        options_storage.args.disable_gzip_output = False
    if options_storage.args.disable_rr is None:
        options_storage.args.disable_rr = False
    if options_storage.args.careful is None:
        options_storage.args.careful = False
    if options_storage.args.mismatch_corrector is None:
        options_storage.args.mismatch_corrector = False
    if options_storage.args.checkpoints is None:
        options_storage.args.checkpoints = "none"
    if options_storage.args.developer_mode is None:
        options_storage.args.developer_mode = False
    if options_storage.args.time_tracer is None:
        options_storage.args.time_tracer = False        
    if options_storage.args.qvoffset == "auto":
        options_storage.args.qvoffset = None
    if options_storage.args.cov_cutoff is None:
        options_storage.args.cov_cutoff = "off"
    if options_storage.args.tmp_dir is None:
        options_storage.args.tmp_dir = os.path.join(options_storage.args.output_dir, options_storage.TMP_DIR)
    if options_storage.args.large_genome is None:
        options_storage.args.large_genome = False
    if options_storage.args.truseq_mode is None:
        options_storage.args.truseq_mode = False
    if options_storage.args.save_gp is None:
        options_storage.args.save_gp = False
    if options_storage.args.only_assembler is None:
        options_storage.args.only_assembler = False
    if options_storage.args.only_error_correction is None:
        options_storage.args.only_error_correction = False


def save_restart_options():
    options_storage.restart = argparse.Namespace(**vars(options_storage.args))
    options_storage.restart.continue_mode = None
    options_storage.restart.restart_from = None
    options_storage.restart.output_dir = None


def load_restart_options():
    if "k_mers" in options_storage.restart and options_storage.restart.k_mers:
        options_storage.original_k_mers = options_storage.args.k_mers
        if options_storage.restart.k_mers == "auto":
            options_storage.args.k_mers = None  # set by default
        else:
            options_storage.args.k_mers = options_storage.restart.k_mers
        options_storage.restart.k_mers = None

    for option in options_storage.restart.__dict__:
        if options_storage.restart.__dict__[option] is not None:
            options_storage.args.__dict__[option] = options_storage.restart.__dict__[option]


def enable_truseq_mode():
    options_storage.K_MERS_SHORT = [21, 33, 45, 55]
    options_storage.K_MERS_150 = [21, 33, 45, 55, 77]
    options_storage.K_MERS_250 = [21, 33, 45, 55, 77, 99, 127]
    options_storage.args.truseq_mode = True
    options_storage.correct_scaffolds = True
    options_storage.run_truseq_postprocessing = True
    options_storage.args.only_assembler = True


def will_rerun(options):
    for opt, arg in options:
        if opt == "--continue" or opt.startswith(
                "--restart-from"):  # checks both --restart-from k33 and --restart-from=k33
            return True
    return False


def is_first_run():
    continue_parser = argparse.ArgumentParser(add_help=False)
    continue_parser.add_argument("--continue", dest="continue_mode", action="store_true")
    continue_parser.add_argument("--restart-from", dest="restart_from", default=None, type=restart_from, action="store")
    nargs, unknown_args = continue_parser.parse_known_args()
    return not (nargs.continue_mode or nargs.restart_from is not None)


def get_output_dir_from_args():
    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("-o", type=str, dest="output_dir", action=StoreUniqueAction)
    nargs, unknown_args = output_parser.parse_known_args()
    if nargs.output_dir is None:
        return None
    return abspath(expanduser(nargs.output_dir))
