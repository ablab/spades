#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
from site import addsitedir
from distutils import dir_util
from os.path import abspath, expanduser
import sys
import getopt
import logging
import platform
import errno

import spades_init
spades_init.init()
spades_home = spades_init.spades_home
bin_home = spades_init.bin_home
python_modules_home = spades_init.python_modules_home
ext_python_modules_home = spades_init.ext_python_modules_home
spades_version = spades_init.spades_version

import support
support.check_python_version()

from process_cfg import merge_configs, empty_config, load_config_from_file
import hammer_logic
import spades_logic
import options_storage
addsitedir(ext_python_modules_home)
if sys.version.startswith('2.'):
    import pyyaml2 as pyyaml
elif sys.version.startswith('3.'):
    import pyyaml3 as pyyaml

import moleculo_postprocessing
import alignment


def print_used_values(cfg, log):
    def print_value(cfg, section, param, pretty_param="", margin="  "):
        if not pretty_param:
            pretty_param = param.capitalize().replace('_', ' ')
        line = margin + pretty_param
        if param in cfg[section].__dict__:
            line += ": " + str(cfg[section].__dict__[param])
        else:
            if param.find("offset") != -1:
                line += " will be auto-detected"
        log.info(line)

    log.info("")

    # system info
    log.info("System information:")
    try:
        log.info("  SPAdes version: " + str(spades_version).strip())
        log.info("  Python version: " + ".".join(map(str, sys.version_info[0:3])))
        # for more details: '[' + str(sys.version_info) + ']'
        log.info("  OS: " + platform.platform())
        # for more details: '[' + str(platform.uname()) + ']'
    except Exception:
        log.info("  Problem occurred when getting system information")
    log.info("")

    # main
    print_value(cfg, "common", "output_dir", "", "")
    if ("error_correction" in cfg) and (not "assembly" in cfg):
        log.info("Mode: ONLY read error correction (without assembling)")
    elif (not "error_correction" in cfg) and ("assembly" in cfg):
        log.info("Mode: ONLY assembling (without read error correction)")
    else:
        log.info("Mode: read error correction and assembling")
    if ("common" in cfg) and ("developer_mode" in cfg["common"].__dict__):
        if cfg["common"].developer_mode:
            log.info("Debug mode is turned ON")
        else:
            log.info("Debug mode is turned OFF")
    log.info("")

    # dataset
    if "dataset" in cfg:
        log.info("Dataset parameters:")

        if options_storage.iontorrent:
            log.info("  IonTorrent data")

        if options_storage.meta:
            log.info("  Metagenomic mode")
        elif options_storage.large_genome:
            log.info("  Large genome mode")
        elif options_storage.truseq_mode:
            log.info("  Illumina TruSeq mode")
        elif options_storage.rna:
            log.info("  RNA-seq mode")
        elif options_storage.single_cell:
            log.info("  Single-cell mode")
        else:
            log.info("  Multi-cell mode (you should set '--sc' flag if input data"\
                     " was obtained with MDA (single-cell) technology"\
                     " or --meta flag if processing metagenomic dataset)")

        log.info("  Reads:")
        dataset_data = pyyaml.load(open(cfg["dataset"].yaml_filename, 'r'))
        dataset_data = support.relative2abs_paths(dataset_data, os.path.dirname(cfg["dataset"].yaml_filename))
        support.pretty_print_reads(dataset_data, log)

    # error correction
    if "error_correction" in cfg:
        log.info("Read error correction parameters:")
        print_value(cfg, "error_correction", "max_iterations", "Iterations")
        print_value(cfg, "error_correction", "qvoffset", "PHRED offset")

        if cfg["error_correction"].gzip_output:
            log.info("  Corrected reads will be compressed (with gzip)")
        else:
            log.info("  Corrected reads will NOT be compressed (with gzip)")

    # assembly
    if "assembly" in cfg:
        log.info("Assembly parameters:")
        if options_storage.auto_K_allowed():
            log.info("  k: automatic selection based on read length")
        else:
            print_value(cfg, "assembly", "iterative_K", "k")
        if options_storage.plasmid:
            log.info("  Plasmid mode is turned ON")
        if cfg["assembly"].disable_rr:
            log.info("  Repeat resolution is DISABLED")
        else:
            log.info("  Repeat resolution is enabled")
        if options_storage.careful:
            log.info("  Mismatch careful mode is turned ON")
        else:
            log.info("  Mismatch careful mode is turned OFF")
        if "mismatch_corrector" in cfg:
            log.info("  MismatchCorrector will be used")
        else:
            log.info("  MismatchCorrector will be SKIPPED")
        if cfg["assembly"].cov_cutoff == 'off':
            log.info("  Coverage cutoff is turned OFF")
        elif cfg["assembly"].cov_cutoff == 'auto':
            log.info("  Coverage cutoff is turned ON and threshold will be auto-detected")
        else:
            log.info("  Coverage cutoff is turned ON and threshold is " + str(cfg["assembly"].cov_cutoff))

    log.info("Other parameters:")
    print_value(cfg, "common", "tmp_dir", "Dir for temp files")
    print_value(cfg, "common", "max_threads", "Threads")
    print_value(cfg, "common", "max_memory", "Memory limit (in Gb)", "  ")
    log.info("")


def fill_cfg(options_to_parse, log, secondary_filling=False):
    skip_output_dir=secondary_filling
    skip_stop_after = secondary_filling
    load_processed_dataset=secondary_filling

    try:
        options, not_options = getopt.gnu_getopt(options_to_parse, options_storage.short_options, options_storage.long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        sys.stderr.flush()
        show_usage(1)

    if not options:
        show_usage(1)

    if len(not_options) > 1:
        for opt, arg in options:
            if opt == "-k" and arg.strip().endswith(','):
                support.error("Do not put spaces after commas in the list of k-mers sizes! Correct example: -k 21,33,55", log)
        support.error("Please specify option (e.g. -1, -2, -s, etc) for the following paths: " + ", ".join(not_options[1:]) + "\n", log)

    # all parameters are stored here
    cfg = dict()
    # dataset is stored here. We are prepared for up to MAX_LIBS_NUMBER for each type of short-reads libs
    dataset_data = [{} for i in range(options_storage.MAX_LIBS_NUMBER *
                                      len(options_storage.SHORT_READS_TYPES.keys()) +
                                      len(options_storage.LONG_READS_TYPES))]  # "[{}]*num" doesn't work here!

    # auto detecting SPAdes mode (rna, meta, etc) if it is not a rerun (--continue or --restart-from)
    if secondary_filling or not options_storage.will_rerun(options):
        mode = options_storage.get_mode()
        if mode is not None:
            options.append(('--' + mode, ''))

    # for parsing options from "previous run command"
    options_storage.continue_mode = False
    options_storage.k_mers = None
    for opt, arg in options:
        if opt == '-o':
            if not skip_output_dir:
                if options_storage.output_dir is not None:
                    support.error('-o option was specified at least twice')
                options_storage.output_dir = abspath(expanduser(arg))
                options_storage.dict_of_rel2abs[arg] = options_storage.output_dir
                support.check_path_is_ascii(options_storage.output_dir, 'output directory')
        elif opt == "--tmp-dir":
            options_storage.tmp_dir = abspath(expanduser(arg))
            options_storage.dict_of_rel2abs[arg] = options_storage.tmp_dir
            support.check_path_is_ascii(options_storage.tmp_dir, 'directory for temporary files')
        elif opt == "--configs-dir":
            options_storage.configs_dir = support.check_dir_existence(arg)
        elif opt == "--reference":
            options_storage.reference = support.check_file_existence(arg, 'reference', log)
            options_storage.developer_mode = True
        elif opt == "--series-analysis":
            options_storage.series_analysis = support.check_file_existence(arg, 'series-analysis', log)
        elif opt == "--dataset":
            options_storage.dataset_yaml_filename = support.check_file_existence(arg, 'dataset', log)

        elif opt in options_storage.reads_options:
            support.add_to_dataset(opt, arg, dataset_data)

        elif opt == '-k':
            if arg == 'auto':
                options_storage.k_mers = arg
            else:
                options_storage.k_mers = list(map(int, arg.split(",")))
                for k in options_storage.k_mers:
                    if k < options_storage.MIN_K or k > options_storage.MAX_K:
                        support.error('wrong k value ' + str(k) + ': all k values should be between %d and %d' %
                                                                  (options_storage.MIN_K, options_storage.MAX_K), log)
                    if k % 2 == 0:
                        support.error('wrong k value ' + str(k) + ': all k values should be odd', log)

        elif opt == "--sc":
            options_storage.single_cell = True
        elif opt == "--meta":
            options_storage.meta = True
        elif opt == "--large-genome":
            options_storage.large_genome = True
        elif opt == "--plasmid":
            options_storage.plasmid = True

        elif opt == "--rna":
            options_storage.rna = True
        elif opt.startswith("--ss-"):  # strand specificity, RNA-Seq only
            if opt == "--ss-rf":
                options_storage.strand_specific = True
            elif opt == "--ss-fr":
                options_storage.strand_specific = False
        elif opt == "--fast":  # fast run, RNA-Seq only
            options_storage.fast = True
        elif opt == "--fast:false":
            options_storage.fast = False

        elif opt == "--iontorrent":
            options_storage.iontorrent = True
        elif opt == "--disable-gzip-output":
            options_storage.disable_gzip_output = True
        elif opt == "--disable-gzip-output:false":
            options_storage.disable_gzip_output = False
        elif opt == "--disable-rr":
            options_storage.disable_rr = True
        elif opt == "--disable-rr:false":
            options_storage.disable_rr = False
        
        elif opt == "--only-error-correction":
            if options_storage.only_assembler:
                support.error('you cannot specify --only-error-correction and --only-assembler simultaneously')
            options_storage.only_error_correction = True
        elif opt == "--only-assembler":
            if options_storage.only_error_correction:
                support.error('you cannot specify --only-error-correction and --only-assembler simultaneously')
            options_storage.only_assembler = True

        elif opt == "--read-buffer-size":
            options_storage.read_buffer_size = int(arg)
        elif opt == "--bh-heap-check":
            options_storage.bh_heap_check = arg
        elif opt == "--spades-heap-check":
            options_storage.spades_heap_check = arg

        elif opt == "--continue":
            options_storage.continue_mode = True
        elif opt == "--restart-from":
            if arg not in ['ec', 'as', 'mc', 'scc', 'tpp'] and not arg.startswith('k'):
                support.error("wrong value for --restart-from option: " + arg +
                              " (should be 'ec', 'as', 'k<int>', or 'mc'", log)
            options_storage.continue_mode = True
            options_storage.restart_from = arg
        elif opt == "--stop-after":
            if not skip_stop_after:
                if arg not in ['ec', 'as', 'mc', 'scc', 'tpp'] and not arg.startswith('k'):
                    support.error("wrong value for --stop-after option: " + arg +
                                  " (should be 'ec', 'as', 'k<int>', or 'mc'", log)
                options_storage.stop_after = arg

        elif opt == '-t' or opt == "--threads":
            options_storage.threads = int(arg)
        elif opt == '-m' or opt == "--memory":
            options_storage.memory = int(arg)
        elif opt == "--phred-offset":
            if arg == 'auto':
                options_storage.qvoffset = arg
            elif arg in ['33', '64']:
                options_storage.qvoffset = int(arg)
            else:
                support.error('wrong PHRED quality offset value: ' + arg +
                              ' (should be either 33, 64, or \'auto\')', log)
        elif opt == "--save-gp":
            options_storage.save_gp = True
        elif opt == "--cov-cutoff":
            if arg == 'auto' or arg == 'off':
                options_storage.cov_cutoff = arg
            elif support.is_float(arg) and float(arg) > 0.0:
                options_storage.cov_cutoff = float(arg)
            else:
                support.error('wrong value for --cov-cutoff option: ' + arg +
                              ' (should be a positive float number, or \'auto\', or \'off\')', log)
        elif opt == "--hidden-cov-cutoff":
            if support.is_float(arg) and float(arg) > 0.0:
                options_storage.lcer_cutoff = float(arg)
            else:
                support.error('wrong value for --hidden-cov-cutoff option: ' + arg +
                              ' (should be a positive float number)', log)
        elif opt == '-i' or opt == "--iterations":
            options_storage.iterations = int(arg)

        elif opt == "--debug":
            options_storage.developer_mode = True
        elif opt == "--debug:false":
            options_storage.developer_mode = False

        #corrector
        elif opt == "--mismatch-correction":
            options_storage.mismatch_corrector = True
        elif opt == "--mismatch-correction:false":
            options_storage.mismatch_corrector = False

        elif opt == "--careful":
            options_storage.mismatch_corrector = True
            options_storage.careful = True
        elif opt == "--careful:false":
            options_storage.mismatch_corrector = False
            options_storage.careful = False

        elif opt == '-v' or opt == "--version":
            show_version()
        elif opt == '-h' or opt == "--help":
            show_usage(0)
        elif opt == "--help-hidden":
            show_usage(0, show_hidden=True)

        elif opt == "--test":
            options_storage.set_test_options()            
            #break
        elif opt == "--diploid":
            options_storage.diploid_mode = True
        elif opt == "--truseq":
            options_storage.enable_truseq_mode()
        else:
            raise ValueError

    if options_storage.test_mode:
        if options_storage.plasmid:
            support.add_to_dataset('-1', os.path.join(spades_home, "test_dataset_plasmid/pl1.fq.gz"), dataset_data)
            support.add_to_dataset('-2', os.path.join(spades_home, "test_dataset_plasmid/pl2.fq.gz"), dataset_data)
        else:
            support.add_to_dataset('-1', os.path.join(spades_home, "test_dataset/ecoli_1K_1.fq.gz"), dataset_data)
            support.add_to_dataset('-2', os.path.join(spades_home, "test_dataset/ecoli_1K_2.fq.gz"), dataset_data)

    if not options_storage.output_dir:
        support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).", log)
    if not os.path.isdir(options_storage.output_dir):
        if options_storage.continue_mode:
            support.error("the output_dir should exist for --continue and for --restart-from!", log)
        os.makedirs(options_storage.output_dir)
    if options_storage.restart_from:
        if options_storage.continue_mode:  # saving parameters specified with --restart-from
            if not support.dataset_is_empty(dataset_data):
                support.error("you cannot specify reads with --restart-from option!", log)
            options_storage.save_restart_options(log)
        else:  # overriding previous run parameters
            options_storage.load_restart_options()
    if options_storage.meta:
        if options_storage.careful or options_storage.mismatch_corrector or options_storage.cov_cutoff != "off":
            support.error("you cannot specify --careful, --mismatch-correction or --cov-cutoff in metagenomic mode!", log)
    if options_storage.rna:
        if options_storage.careful:
            support.error("you cannot specify --careful in RNA-Seq mode!", log)
        if options_storage.k_mers and options_storage.k_mers != 'auto' and len(options_storage.k_mers) > 1:
            support.error("you cannot specify multiple k-mer sizes in RNA-Seq mode!", log)
    if [options_storage.meta, options_storage.large_genome, options_storage.truseq_mode,
       options_storage.rna, options_storage.plasmid, options_storage.single_cell].count(True) > 1:
        support.error("you cannot simultaneously use more than one mode out of "
                      "Metagenomic, Large genome, Illumina TruSeq, RNA-Seq, Plasmid, and Single-cell!", log)
    if options_storage.continue_mode:
        return None, None

    existing_dataset_data = None
    processed_dataset_fpath = os.path.join(options_storage.output_dir, "input_dataset.yaml")
    if load_processed_dataset:
        if os.path.isfile(processed_dataset_fpath):
            try:
                existing_dataset_data = pyyaml.load(open(processed_dataset_fpath, 'r'))
            except pyyaml.YAMLError:
                existing_dataset_data = None
    if existing_dataset_data is not None:
        dataset_data = existing_dataset_data
    else:
        if options_storage.dataset_yaml_filename:
            try:
                dataset_data = pyyaml.load(open(options_storage.dataset_yaml_filename, 'r'))
            except pyyaml.YAMLError:
                _, exc, _ = sys.exc_info()
                support.error('exception caught while parsing YAML file (' + options_storage.dataset_yaml_filename + '):\n' + str(exc))
            dataset_data = support.relative2abs_paths(dataset_data, os.path.dirname(options_storage.dataset_yaml_filename))
        else:
            dataset_data = support.correct_dataset(dataset_data)
            dataset_data = support.relative2abs_paths(dataset_data, os.getcwd())
    options_storage.dataset_yaml_filename = processed_dataset_fpath

    support.check_dataset_reads(dataset_data, options_storage.only_assembler, log)
    if not support.get_lib_ids_by_type(dataset_data, spades_logic.READS_TYPES_USED_IN_CONSTRUCTION):
        support.error('you should specify at least one unpaired, paired-end, or high-quality mate-pairs library!')
    if options_storage.rna:
        if len(dataset_data) != len(support.get_lib_ids_by_type(dataset_data, spades_logic.READS_TYPES_USED_IN_RNA_SEQ)):
            support.error('you cannot specify any data types except ' +
                          ', '.join(spades_logic.READS_TYPES_USED_IN_RNA_SEQ) + ' in RNA-Seq mode!')
        #if len(support.get_lib_ids_by_type(dataset_data, 'paired-end')) > 1:
        #    support.error('you cannot specify more than one paired-end library in RNA-Seq mode!')

    if existing_dataset_data is None:
        pyyaml.dump(dataset_data, open(options_storage.dataset_yaml_filename, 'w'),
                    default_flow_style=False, default_style='"', width=float("inf"))

    options_storage.set_default_values()
    ### FILLING cfg
    cfg["common"] = empty_config()
    cfg["dataset"] = empty_config()
    if not options_storage.only_assembler:
        cfg["error_correction"] = empty_config()
    if not options_storage.only_error_correction:
        cfg["assembly"] = empty_config()

    # common
    cfg["common"].__dict__["output_dir"] = options_storage.output_dir
    cfg["common"].__dict__["tmp_dir"] = options_storage.tmp_dir
    cfg["common"].__dict__["max_threads"] = options_storage.threads
    cfg["common"].__dict__["max_memory"] = options_storage.memory
    cfg["common"].__dict__["developer_mode"] = options_storage.developer_mode
    if options_storage.series_analysis:
        cfg["common"].__dict__["series_analysis"] = options_storage.series_analysis

    # dataset section
    cfg["dataset"].__dict__["yaml_filename"] = options_storage.dataset_yaml_filename
    if options_storage.developer_mode and options_storage.reference:
        cfg["dataset"].__dict__["reference"] = options_storage.reference

    # error correction
    if (not options_storage.only_assembler) and (options_storage.iterations > 0):
        cfg["error_correction"].__dict__["output_dir"] = os.path.join(cfg["common"].output_dir, "corrected")
        cfg["error_correction"].__dict__["max_iterations"] = options_storage.iterations
        cfg["error_correction"].__dict__["gzip_output"] = not options_storage.disable_gzip_output
        if options_storage.qvoffset:
            cfg["error_correction"].__dict__["qvoffset"] = options_storage.qvoffset
        if options_storage.bh_heap_check:
            cfg["error_correction"].__dict__["heap_check"] = options_storage.bh_heap_check
        cfg["error_correction"].__dict__["iontorrent"] = options_storage.iontorrent
        if options_storage.meta or options_storage.large_genome:
            cfg["error_correction"].__dict__["count_filter_singletons"] = 1
        if options_storage.read_buffer_size:
            cfg["error_correction"].__dict__["read_buffer_size"] = options_storage.read_buffer_size

    # assembly
    if not options_storage.only_error_correction:
        if options_storage.k_mers == 'auto' and options_storage.restart_from is None:
            options_storage.k_mers = None
        if options_storage.k_mers:
            cfg["assembly"].__dict__["iterative_K"] = options_storage.k_mers
        elif options_storage.rna:
            k_value = options_storage.K_MERS_RNA
            if not options_storage.iontorrent:
                k_value = support.get_reads_length(dataset_data, log) / 2 - 1
                if k_value % 2 == 0:
                    k_value -= 1
            cfg["assembly"].__dict__["iterative_K"] = k_value
        else:
            cfg["assembly"].__dict__["iterative_K"] = options_storage.K_MERS_SHORT
        cfg["assembly"].__dict__["disable_rr"] = options_storage.disable_rr
        cfg["assembly"].__dict__["diploid_mode"] = options_storage.diploid_mode
        cfg["assembly"].__dict__["cov_cutoff"] = options_storage.cov_cutoff
        cfg["assembly"].__dict__["lcer_cutoff"] = options_storage.lcer_cutoff
        cfg["assembly"].__dict__["save_gp"] = options_storage.save_gp
        if options_storage.spades_heap_check:
            cfg["assembly"].__dict__["heap_check"] = options_storage.spades_heap_check
        if options_storage.read_buffer_size:
            cfg["assembly"].__dict__["read_buffer_size"] = options_storage.read_buffer_size
        cfg["assembly"].__dict__["correct_scaffolds"] = options_storage.correct_scaffolds

    #corrector can work only if contigs exist (not only error correction)
    if (not options_storage.only_error_correction) and options_storage.mismatch_corrector:
        cfg["mismatch_corrector"] = empty_config()
        cfg["mismatch_corrector"].__dict__["skip-masked"] = None
        cfg["mismatch_corrector"].__dict__["bwa"] = os.path.join(bin_home, "bwa-spades")
        cfg["mismatch_corrector"].__dict__["threads"] = options_storage.threads
        cfg["mismatch_corrector"].__dict__["output-dir"] = options_storage.output_dir
    cfg["run_truseq_postprocessing"] = options_storage.run_truseq_postprocessing
    return cfg, dataset_data

def check_cfg_for_partial_run(cfg, type='restart-from'):  # restart-from ot stop-after
    if type == 'restart-from':
        check_point = options_storage.restart_from
        action = 'restart from'
        verb = 'was'
    elif type == 'stop-after':
        check_point = options_storage.stop_after
        action = 'stop after'
        verb = 'is'
    else:
        return

    if check_point == 'ec' and ("error_correction" not in cfg):
        support.error("failed to " + action + " 'read error correction' ('" + check_point + "') because this stage " + verb + " not specified!")
    if check_point == 'mc' and ("mismatch_corrector" not in cfg):
        support.error("failed to " + action + " 'mismatch correction' ('" + check_point + "') because this stage " + verb + " not specified!")
    if check_point == 'as' or check_point.startswith('k'):
        if "assembly" not in cfg:
            support.error("failed to " + action + " 'assembling' ('" + check_point + "') because this stage " + verb + " not specified!")
        if check_point.startswith('k'):
            correct_k = False
            k_to_check = options_storage.k_mers
            if not k_to_check:
                if options_storage.auto_K_allowed():
                    k_to_check = list(set(options_storage.K_MERS_SHORT + options_storage.K_MERS_150 + options_storage.K_MERS_250))
                else:
                    k_to_check = options_storage.K_MERS_SHORT
            for k in k_to_check:
                if check_point == ("k%d" % k) or check_point.startswith("k%d:" % k):
                    correct_k = True
                    break
            if not correct_k:
                k_str = check_point[1:]
                if k_str.find(":") != -1:
                    k_str = k_str[:k_str.find(":")]
                support.error("failed to " + action + " K=%s because this K " % k_str + verb + " not specified!")


def get_options_from_params(params_filename, running_script):
    cmd_line = None
    options = None
    if not os.path.isfile(params_filename):
        return cmd_line, options, "failed to parse command line of the previous run (%s not found)!" % params_filename
    params = open(params_filename, 'r')
    cmd_line = params.readline().strip()
    spades_prev_version = None
    for line in params:
        if line.find('SPAdes version:') != -1:
            spades_prev_version = line.split('SPAdes version:')[1]
            break
    params.close()
    if spades_prev_version is None:
        return cmd_line, options, "failed to parse SPAdes version of the previous run!"
    if spades_prev_version.strip() != spades_version.strip():
        return cmd_line, options, "SPAdes version of the previous run (%s) is not equal " \
                                  "to the current version of SPAdes (%s)!" \
                                  % (spades_prev_version.strip(), spades_version.strip())
    if 'Command line: ' not in cmd_line or '\t' not in cmd_line:
        return cmd_line, options, "failed to parse executable script of the previous run!"
    options = cmd_line.split('\t')[1:]
    prev_running_script = cmd_line.split('\t')[0][len('Command line: '):]
    # we cannot restart/continue spades.py run with metaspades.py/rnaspades.py/etc and vice versa
    if os.path.basename(prev_running_script) != os.path.basename(running_script):
        return cmd_line, options, "executable script of the previous run (%s) is not equal " \
                                  "to the current executable script (%s)!" \
                                  % (os.path.basename(prev_running_script),
                                     os.path.basename(running_script))
    return cmd_line, options, ""


def show_version():
    options_storage.version(spades_version)
    sys.exit(0)


def show_usage(code, show_hidden=False):
    options_storage.usage(spades_version, show_hidden=show_hidden)
    sys.exit(code)


def main(args):
    os.environ["LC_ALL"] = "C"

    if len(args) == 1:
        show_usage(0)

    log = logging.getLogger('spades')
    log.setLevel(logging.DEBUG)

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    support.check_binaries(bin_home, log)

    # parse options and safe all parameters to cfg
    options = args
    cfg, dataset_data = fill_cfg(options, log)

    if options_storage.continue_mode:
        cmd_line, options, err_msg = get_options_from_params(os.path.join(options_storage.output_dir, "params.txt"), args[0])
        if err_msg:
            support.error(err_msg + " Please restart from the beginning or specify another output directory.")
        cfg, dataset_data = fill_cfg(options, log, secondary_filling=True)
        if options_storage.restart_from:
            check_cfg_for_partial_run(cfg, type='restart-from')
        options_storage.continue_mode = True
    if options_storage.stop_after:
        check_cfg_for_partial_run(cfg, type='stop-after')

    log_filename = os.path.join(cfg["common"].output_dir, "spades.log")
    if options_storage.continue_mode:
        log_handler = logging.FileHandler(log_filename, mode='a')
    else:
        log_handler = logging.FileHandler(log_filename, mode='w')
    log.addHandler(log_handler)

    if options_storage.continue_mode:
        log.info("\n======= SPAdes pipeline continued. Log can be found here: " + log_filename + "\n")
        log.info("Restored from " + cmd_line)
        if options_storage.restart_from:
            updated_params = ""
            skip_next = False
            for v in args[1:]:
                if v == '-o' or v == '--restart-from':
                    skip_next = True
                    continue
                if skip_next or v.startswith('--restart-from='):  # you can specify '--restart-from=k33' but not '-o=out_dir'
                    skip_next = False
                    continue
                updated_params += "\t" + v
            updated_params = updated_params.strip()
            log.info("with updated parameters: " + updated_params)
            cmd_line += "\t" + updated_params
        log.info("")

    params_filename = os.path.join(cfg["common"].output_dir, "params.txt")
    params_handler = logging.FileHandler(params_filename, mode='w')
    log.addHandler(params_handler)

    if options_storage.continue_mode:
        log.info(cmd_line)
    else:
        command = "Command line: "
        for v in args:
            # substituting relative paths with absolute ones (read paths, output dir path, etc)
            v, prefix = support.get_option_prefix(v)
            if v in options_storage.dict_of_rel2abs.keys():
                v = options_storage.dict_of_rel2abs[v]
            if prefix:
                command += prefix + ":"
            command += v + "\t"
        log.info(command)

    # special case
#    if "mismatch_corrector" in cfg and not support.get_lib_ids_by_type(dataset_data, 'paired-end'):
#        support.warning('cannot perform mismatch correction without at least one paired-end library! Skipping this step.', log)
#        del cfg["mismatch_corrector"]

    print_used_values(cfg, log)
    log.removeHandler(params_handler)

    support.check_single_reads_in_options(options, log)

    if not options_storage.continue_mode:
        log.info("\n======= SPAdes pipeline started. Log can be found here: " + log_filename + "\n")

    # splitting interlaced reads and processing Ns in additional contigs if needed
    if support.dataset_has_interlaced_reads(dataset_data) or support.dataset_has_additional_contigs(dataset_data)\
            or support.dataset_has_nxmate_reads(dataset_data):
        dir_for_split_reads = os.path.join(options_storage.output_dir, 'split_input')
        if support.dataset_has_interlaced_reads(dataset_data) or support.dataset_has_nxmate_reads(dataset_data):
            if not os.path.isdir(dir_for_split_reads):
                os.makedirs(dir_for_split_reads)
            if support.dataset_has_interlaced_reads(dataset_data):
                dataset_data = support.split_interlaced_reads(dataset_data, dir_for_split_reads, log)
            if support.dataset_has_nxmate_reads(dataset_data):
                dataset_data = support.process_nxmate_reads(dataset_data, dir_for_split_reads, log)
        if support.dataset_has_additional_contigs(dataset_data):
            dataset_data = support.process_Ns_in_additional_contigs(dataset_data, dir_for_split_reads, log)
        options_storage.dataset_yaml_filename = os.path.join(options_storage.output_dir, "input_dataset.yaml")
        pyyaml.dump(dataset_data, open(options_storage.dataset_yaml_filename, 'w'),
                    default_flow_style=False, default_style='"', width=float("inf"))
        cfg["dataset"].yaml_filename = options_storage.dataset_yaml_filename

    try:
        # copying configs before all computations (to prevent its changing at run time)
        tmp_configs_dir = os.path.join(cfg["common"].output_dir, "configs")
        if os.path.isdir(tmp_configs_dir) and not options_storage.continue_mode:
            shutil.rmtree(tmp_configs_dir)
        if not os.path.isdir(tmp_configs_dir):
            if options_storage.configs_dir:
                dir_util.copy_tree(options_storage.configs_dir, tmp_configs_dir, preserve_times=False, preserve_mode=False)
            else:
                dir_util.copy_tree(os.path.join(spades_home, "configs"), tmp_configs_dir, preserve_times=False, preserve_mode=False)

        corrected_dataset_yaml_filename = ''
        if "error_correction" in cfg:
            STAGE_NAME = "Read error correction"
            bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])
            corrected_dataset_yaml_filename = os.path.join(bh_cfg.output_dir, "corrected.yaml")
            ec_is_needed = True
            only_compressing_is_needed = False
            if os.path.isfile(corrected_dataset_yaml_filename) and options_storage.continue_mode \
                    and not options_storage.restart_from == "ec":
                if not bh_cfg.gzip_output or \
                        support.dataset_has_gzipped_reads(pyyaml.load(open(corrected_dataset_yaml_filename, 'r'))):
                    log.info("\n===== Skipping %s (already processed). \n" % STAGE_NAME)
                    ec_is_needed = False
                else:
                    only_compressing_is_needed = True
            if ec_is_needed:
                if not only_compressing_is_needed:
                    support.continue_from_here(log)

                    if "HEAPCHECK" in os.environ:
                        del os.environ["HEAPCHECK"]
                    if "heap_check" in bh_cfg.__dict__:
                        os.environ["HEAPCHECK"] = bh_cfg.heap_check

                    if os.path.exists(bh_cfg.output_dir):
                        shutil.rmtree(bh_cfg.output_dir)
                    os.makedirs(bh_cfg.output_dir)

                bh_cfg.__dict__["dataset_yaml_filename"] = cfg["dataset"].yaml_filename
                log.info("\n===== %s started. \n" % STAGE_NAME)

                hammer_logic.run_hammer(corrected_dataset_yaml_filename, tmp_configs_dir, bin_home, bh_cfg, dataset_data,
                    ext_python_modules_home, only_compressing_is_needed, log)
                log.info("\n===== %s finished. \n" % STAGE_NAME)
            if options_storage.stop_after == 'ec':
                support.finish_here(log)

        result_contigs_filename = os.path.join(cfg["common"].output_dir, options_storage.contigs_name)
        result_scaffolds_filename = os.path.join(cfg["common"].output_dir, options_storage.scaffolds_name)
        result_assembly_graph_filename = os.path.join(cfg["common"].output_dir, options_storage.assembly_graph_name)
        result_assembly_graph_filename_gfa = os.path.join(cfg["common"].output_dir, options_storage.assembly_graph_name_gfa)
        result_contigs_paths_filename = os.path.join(cfg["common"].output_dir, options_storage.contigs_paths)
        result_scaffolds_paths_filename = os.path.join(cfg["common"].output_dir, options_storage.scaffolds_paths)
        result_transcripts_filename = os.path.join(cfg["common"].output_dir, options_storage.transcripts_name)
        result_transcripts_paths_filename = os.path.join(cfg["common"].output_dir, options_storage.transcripts_paths)
        truseq_long_reads_file_base = os.path.join(cfg["common"].output_dir, "truseq_long_reads")
        truseq_long_reads_file = truseq_long_reads_file_base + ".fasta"
        misc_dir = os.path.join(cfg["common"].output_dir, "misc")
        ### if mismatch correction is enabled then result contigs are copied to misc directory
        assembled_contigs_filename = os.path.join(misc_dir, "assembled_contigs.fasta")
        assembled_scaffolds_filename = os.path.join(misc_dir, "assembled_scaffolds.fasta")
        if "assembly" in cfg and not options_storage.run_completed:
            STAGE_NAME = "Assembling"
            spades_cfg = merge_configs(cfg["assembly"], cfg["common"])
            spades_cfg.__dict__["result_contigs"] = result_contigs_filename
            spades_cfg.__dict__["result_scaffolds"] = result_scaffolds_filename
            spades_cfg.__dict__["result_graph"] = result_assembly_graph_filename
            spades_cfg.__dict__["result_graph_gfa"] = result_assembly_graph_filename_gfa
            spades_cfg.__dict__["result_contigs_paths"] = result_contigs_paths_filename
            spades_cfg.__dict__["result_scaffolds_paths"] = result_scaffolds_paths_filename
            spades_cfg.__dict__["result_transcripts"] = result_transcripts_filename
            spades_cfg.__dict__["result_transcripts_paths"] = result_transcripts_paths_filename

            if options_storage.continue_mode and (os.path.isfile(spades_cfg.result_contigs)
                                                  or ("mismatch_corrector" in cfg and
                                                      os.path.isfile(assembled_contigs_filename))
                                                  or (options_storage.truseq_mode and os.path.isfile(assembled_scaffolds_filename)))\
                and not options_storage.restart_from == 'as' \
                and not options_storage.restart_from == 'scc' \
                and not (options_storage.restart_from and options_storage.restart_from.startswith('k')):

                log.info("\n===== Skipping %s (already processed). \n" % STAGE_NAME)
                # calculating latest_dir for the next stages
                latest_dir = support.get_latest_dir(os.path.join(spades_cfg.output_dir, "K*"))
                if not latest_dir:
                    support.error("failed to continue the previous run! Please restart from previous stages or from the beginning.", log)
            else:
                old_result_files = [result_contigs_filename, result_scaffolds_filename,
                                    assembled_contigs_filename, assembled_scaffolds_filename]
                for old_result_file in old_result_files:
                    if os.path.isfile(old_result_file):
                        os.remove(old_result_file)

                if options_storage.restart_from == 'as':
                    support.continue_from_here(log)

                if os.path.isfile(corrected_dataset_yaml_filename):
                    dataset_data = pyyaml.load(open(corrected_dataset_yaml_filename, 'r'))
                    dataset_data = support.relative2abs_paths(dataset_data, os.path.dirname(corrected_dataset_yaml_filename))
                if spades_cfg.disable_rr:
                    spades_cfg.__dict__["rr_enable"] = False
                else:
                    spades_cfg.__dict__["rr_enable"] = True

                if "HEAPCHECK" in os.environ:
                    del os.environ["HEAPCHECK"]
                if "heap_check" in spades_cfg.__dict__:
                    os.environ["HEAPCHECK"] = spades_cfg.heap_check

                log.info("\n===== %s started.\n" % STAGE_NAME)

                # creating dataset
                dataset_filename = os.path.join(spades_cfg.output_dir, "dataset.info")
                if not os.path.isfile(dataset_filename) or not options_storage.continue_mode:
                    dataset_file = open(dataset_filename, 'w')
                    import process_cfg
                    if os.path.isfile(corrected_dataset_yaml_filename):
                        dataset_file.write("reads" + '\t' + process_cfg.process_spaces(corrected_dataset_yaml_filename) + '\n')
                    else:
                        dataset_file.write("reads" + '\t' + process_cfg.process_spaces(cfg["dataset"].yaml_filename) + '\n')
                    if spades_cfg.developer_mode and "reference" in cfg["dataset"].__dict__:
                        dataset_file.write("reference_genome" + '\t')
                        dataset_file.write(process_cfg.process_spaces(cfg["dataset"].reference) + '\n')
                    dataset_file.close()
                spades_cfg.__dict__["dataset"] = dataset_filename

                used_K = spades_logic.run_spades(tmp_configs_dir, bin_home, spades_cfg, dataset_data, ext_python_modules_home, log)

                if os.path.isdir(misc_dir) and not options_storage.continue_mode:
                    shutil.rmtree(misc_dir)
                if not os.path.isdir(misc_dir):
                    os.makedirs(misc_dir)

                if options_storage.continue_mode and options_storage.restart_from and options_storage.restart_from.startswith('k'):
                    k_str = options_storage.restart_from[1:]
                    if k_str.find(":") != -1:
                        k_str = k_str[:k_str.find(":")]
                    support.error("failed to continue from K=%s because this K was not processed in the original run!" % k_str, log)
                log.info("\n===== %s finished. Used k-mer sizes: %s \n" % (STAGE_NAME, ', '.join(map(str, used_K))))
            if not options_storage.run_completed:
                if options_storage.stop_after == 'as' or options_storage.stop_after == 'scc' or (options_storage.stop_after and options_storage.stop_after.startswith('k')):
                    support.finish_here(log)

            #postprocessing
            if cfg["run_truseq_postprocessing"] and not options_storage.run_completed:
                if options_storage.continue_mode and os.path.isfile(truseq_long_reads_file_base + ".fastq") and not options_storage.restart_from == 'tpp':
                    log.info("\n===== Skipping %s (already processed). \n" % "TruSeq postprocessing")
                else:
                    support.continue_from_here(log)
                    if os.path.isfile(result_scaffolds_filename):
                        shutil.move(result_scaffolds_filename, assembled_scaffolds_filename)
                    reads_library = dataset_data[0]
                    alignment_bin = os.path.join(bin_home, "bwa-spades")
                    alignment_dir = os.path.join(cfg["common"].output_dir, "alignment")
                    sam_files = alignment.align_bwa(alignment_bin, assembled_scaffolds_filename, dataset_data, alignment_dir, log, options_storage.threads)
                    moleculo_postprocessing.moleculo_postprocessing(assembled_scaffolds_filename, truseq_long_reads_file_base, sam_files, log)
                if options_storage.stop_after == 'tpp':
                    support.finish_here(log)

            #corrector
            if "mismatch_corrector" in cfg and not options_storage.run_completed and \
                    (os.path.isfile(result_contigs_filename) or
                    (options_storage.continue_mode and os.path.isfile(assembled_contigs_filename))):
                STAGE_NAME = "Mismatch correction"
                to_correct = dict()
                to_correct["contigs"] = (result_contigs_filename, assembled_contigs_filename)
                if os.path.isfile(result_scaffolds_filename) or (options_storage.continue_mode and
                                                                 os.path.isfile(assembled_scaffolds_filename)):
                    to_correct["scaffolds"] = (result_scaffolds_filename, assembled_scaffolds_filename)

                # moving assembled contigs (scaffolds) to misc dir
                for assembly_type, (old, new) in to_correct.items():
                    if options_storage.continue_mode and os.path.isfile(new):
                        continue
                    if os.path.isfile(old):
                        shutil.move(old, new)

                if options_storage.continue_mode and os.path.isfile(result_contigs_filename) and \
                    (os.path.isfile(result_scaffolds_filename) or not os.path.isfile(assembled_scaffolds_filename)) \
                    and not options_storage.restart_from == 'mc':
                    log.info("\n===== Skipping %s (already processed). \n" % STAGE_NAME)
                else:
                    if options_storage.restart_from == 'mc':
                        support.continue_from_here(log)

                    log.info("\n===== %s started." % STAGE_NAME)
                    # detecting paired-end library with the largest insert size
                    cfg["mismatch_corrector"].__dict__["dataset"] = cfg["dataset"].yaml_filename
                    #TODO: add reads orientation

                    import corrector_logic
                    corrector_cfg = cfg["mismatch_corrector"]
                    # processing contigs and scaffolds (or only contigs)
                    for assembly_type, (corrected, assembled) in to_correct.items():
                        if options_storage.continue_mode and os.path.isfile(corrected):
                            log.info("\n== Skipping processing of " + assembly_type + " (already processed)\n")
                            continue
                        if not os.path.isfile(assembled) or os.path.getsize(assembled) == 0:
                            log.info("\n== Skipping processing of " + assembly_type + " (empty file)\n")
                            continue
                        support.continue_from_here(log)
                        log.info("\n== Processing of " + assembly_type + "\n")

                        tmp_dir_for_corrector = os.path.join(cfg["common"].output_dir, "mismatch_corrector", assembly_type)

                        cfg["mismatch_corrector"].__dict__["output_dir"] = tmp_dir_for_corrector
                        # correcting
                        corr_cfg = merge_configs(cfg["mismatch_corrector"], cfg["common"])
                        
                        result_corrected_filename = os.path.join(tmp_dir_for_corrector, "corrected_contigs.fasta")
                        corrector_logic.run_corrector( tmp_configs_dir, bin_home, corr_cfg,
                        ext_python_modules_home, log, assembled, result_corrected_filename)

                        if os.path.isfile(result_corrected_filename):
                            shutil.copyfile(result_corrected_filename, corrected)
                        tmp_d = os.path.join(tmp_dir_for_corrector, "tmp")
                        if os.path.isdir(tmp_d) and not cfg["common"].developer_mode:
                            shutil.rmtree(tmp_d)
                    log.info("\n===== %s finished.\n" % STAGE_NAME)
                if options_storage.stop_after == 'mc':
                    support.finish_here(log)

        if not cfg["common"].developer_mode and os.path.isdir(tmp_configs_dir):
            shutil.rmtree(tmp_configs_dir)

        if not options_storage.run_completed:
            #log.info("")
            if "error_correction" in cfg and os.path.isdir(os.path.dirname(corrected_dataset_yaml_filename)):
                log.info(" * Corrected reads are in " + support.process_spaces(os.path.dirname(corrected_dataset_yaml_filename) + "/"))
            if "assembly" in cfg and os.path.isfile(result_contigs_filename):
                message = " * Assembled contigs are in " + support.process_spaces(result_contigs_filename)
                log.info(message)
            if options_storage.rna and "assembly" in cfg:
                if os.path.isfile(result_transcripts_filename):
                    message = " * Assembled transcripts are in " + support.process_spaces(result_transcripts_filename)
                    log.info(message)
                if os.path.isfile(result_transcripts_paths_filename):
                    message = " * Paths in the assembly graph corresponding to the transcripts are in " + \
                              support.process_spaces(result_transcripts_paths_filename)
                    log.info(message)
                for filtering_type in options_storage.filtering_types:
                    result_filtered_transcripts_filename = os.path.join(cfg["common"].output_dir,
                                                                        filtering_type + "_filtered_" +
                                                                        options_storage.transcripts_name)
                    if os.path.isfile(result_filtered_transcripts_filename):
                        message = " * " + filtering_type.capitalize() + " filtered transcripts are in " + \
                                  support.process_spaces(result_filtered_transcripts_filename)
                        log.info(message)
            elif "assembly" in cfg:
                if os.path.isfile(result_scaffolds_filename):
                    message = " * Assembled scaffolds are in " + support.process_spaces(result_scaffolds_filename)
                    log.info(message)
                if os.path.isfile(result_assembly_graph_filename):
                    message = " * Assembly graph is in " + support.process_spaces(result_assembly_graph_filename)
                    log.info(message)
                if os.path.isfile(result_assembly_graph_filename_gfa):
                    message = " * Assembly graph in GFA format is in " + support.process_spaces(result_assembly_graph_filename_gfa)
                    log.info(message)
                if os.path.isfile(result_contigs_paths_filename):
                    message = " * Paths in the assembly graph corresponding to the contigs are in " + \
                              support.process_spaces(result_contigs_paths_filename)
                    log.info(message)
                if os.path.isfile(result_scaffolds_paths_filename):
                    message = " * Paths in the assembly graph corresponding to the scaffolds are in " + \
                              support.process_spaces(result_scaffolds_paths_filename)
                    log.info(message)
            #log.info("")

        #breaking scaffolds
        if os.path.isfile(result_scaffolds_filename):
            if not os.path.isdir(misc_dir):
                os.makedirs(misc_dir)
            result_broken_scaffolds = os.path.join(misc_dir, "broken_scaffolds.fasta")
            if not os.path.isfile(result_broken_scaffolds) or not options_storage.continue_mode:
                modified, broken_scaffolds = support.break_scaffolds(result_scaffolds_filename,
                    options_storage.THRESHOLD_FOR_BREAKING_SCAFFOLDS)
                if modified:
                    support.write_fasta(result_broken_scaffolds, broken_scaffolds)
                    #log.info(" * Scaffolds broken by " + str(options_storage.THRESHOLD_FOR_BREAKING_SCAFFOLDS) +
                    # " Ns are in " + result_broken_scaffolds)

        ### printing WARNINGS SUMMARY
        if not support.log_warnings(log):
            log.info("\n======= SPAdes pipeline finished.")  # otherwise it finished WITH WARNINGS

        if options_storage.test_mode:
            if options_storage.truseq_mode:
                if not os.path.isfile(truseq_long_reads_file):
                    support.error("TEST FAILED: %s does not exist!" % truseq_long_reads_file)
            elif options_storage.rna:
                if not os.path.isfile(result_transcripts_filename):
                    support.error("TEST FAILED: %s does not exist!" % result_transcripts_filename)
            else:
                for result_filename in [result_contigs_filename, result_scaffolds_filename]:
                    if os.path.isfile(result_filename):
                        result_fasta = list(support.read_fasta(result_filename))
                        # correctness check: should be one contig of length 1000 bp
                        correct_number = 1
                        if options_storage.plasmid:
                            correct_length = 9667
                        else:
                            correct_length = 1000
                        if not len(result_fasta):
                            support.error("TEST FAILED: %s does not contain contigs!" % result_filename)
                        elif len(result_fasta) > correct_number:
                            support.error("TEST FAILED: %s contains more than %d contig (%d)!" %
                                          (result_filename, correct_number, len(result_fasta)))
                        elif len(result_fasta[0][1]) != correct_length:
                            if len(result_fasta[0][1]) > correct_length:
                                relation = "more"
                            else:
                                relation = "less"
                            support.error("TEST FAILED: %s contains %s than %d bp (%d bp)!" %
                                          (result_filename, relation, correct_length, len(result_fasta[0][1])))
                    else:
                        support.error("TEST FAILED: " + result_filename + " does not exist!")
            log.info("\n========= TEST PASSED CORRECTLY.")


        log.info("\nSPAdes log can be found here: " + log_filename)
        log.info("")
        log.info("Thank you for using SPAdes!")
        log.removeHandler(log_handler)

    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            if exc_type == OSError and exc_value.errno == errno.ENOEXEC: # Exec format error
                support.error("It looks like you are using SPAdes binaries for another platform.\n" +
                              support.get_spades_binaries_info_message())
            else:
                log.exception(exc_value)
                support.error("exception caught: %s" % exc_type, log)
    except BaseException:  # since python 2.5 system-exiting exceptions (e.g. KeyboardInterrupt) are derived from BaseException
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            support.error("exception caught: %s" % exc_type, log)


if __name__ == '__main__':
    main(sys.argv)
