#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
from site import addsitedir
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
        # for more deatils: '[' + str(platform.uname()) + ']'
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

        if cfg["dataset"].single_cell:
            log.info("  Single-cell mode")
        else:
            log.info("  Multi-cell mode (you should set '--sc' flag if input data"\
                     " was obtained with MDA (single-cell) technology")
        if cfg["dataset"].iontorrent:
            log.info("  IonTorrent data")

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
        if cfg["assembly"].careful:
            log.info("  MismatchCorrector will be used")
        else:
            log.info("  MismatchCorrector will be SKIPPED")
        if cfg["assembly"].disable_rr:
            log.info("  Repeat resolution is DISABLED")
        else:
            log.info("  Repeat resolution is enabled")

    log.info("Other parameters:")
    print_value(cfg, "common", "tmp_dir", "Dir for temp files")
    print_value(cfg, "common", "max_threads", "Threads")
    print_value(cfg, "common", "max_memory", "Memory limit (in Gb)", "  ")
    log.info("")


def fill_cfg(options_to_parse, log):
    try:
        options, not_options = getopt.gnu_getopt(options_to_parse, options_storage.short_options, options_storage.long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        sys.stderr.flush()
        options_storage.usage(spades_version)
        sys.exit(1)

    if not options:
        options_storage.usage(spades_version)
        sys.exit(1)

    if len(not_options) > 1:
        support.error("Please specify option (e.g. -1, -2, -s, etc) for the following paths: " + ", ".join(not_options[1:]) + "\n", log)

    # all parameters are stored here
    cfg = dict()
    # dataset is stored here. We are prepared for up to MAX_LIBS_NUMBER paired-end libs and MAX_LIBS_NUMBER mate-pair libs
    dataset_data = [{} for i in range(options_storage.MAX_LIBS_NUMBER * 2)]  # "[{}] * num" doesn't work here!

    options_storage.continue_mode = False
    for opt, arg in options:
        if opt == '-o':
            options_storage.output_dir = os.path.abspath(arg)
        elif opt == "--tmp-dir":
            options_storage.tmp_dir = os.path.abspath(arg)
        elif opt == "--reference":
            options_storage.reference = support.check_file_existence(arg, 'reference', log)
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
            if arg not in ['ec', 'as', 'mc'] and not arg.startswith('k'):
                support.error("wrong value for --restart-from option: " + arg + " (only 'ec', 'as', 'k<int>', 'mc' are available)", log)
            options_storage.continue_mode = True
            options_storage.restart_from = arg

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
                support.error('wrong PHRED quality offset value ' + str(arg) + ': should be either 33 or 64', log)
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

        elif opt == '-h' or opt == "--help":
            options_storage.usage(spades_version)
            sys.exit(0)
        elif opt == "--help-hidden":
            options_storage.usage(spades_version, True)
            sys.exit(0)

        elif opt == "--test":
            options_storage.set_test_options()
            support.add_to_dataset('-1', os.path.join(spades_home, "test_dataset/ecoli_1K_1.fq.gz"), dataset_data)
            support.add_to_dataset('-2', os.path.join(spades_home, "test_dataset/ecoli_1K_2.fq.gz"), dataset_data)
            #break
        elif opt == "--diploid":
            options_storage.diploid_mode = True
        else:
            raise ValueError


    if not options_storage.output_dir:
        support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).", log)
    if not os.path.isdir(options_storage.output_dir):
        if options_storage.continue_mode:
            support.error("the output_dir should exist for --continue and for --restart-from!", log)
        os.makedirs(options_storage.output_dir)
    if options_storage.restart_from:
        if options_storage.continue_mode: # saving parameters specified with --restart-from
            if not support.dataset_is_empty(dataset_data):
                support.error("you cannot specify reads with --restart-from option!", log)
            options_storage.save_restart_options(log)
        else:  # overriding previous run parameters
            options_storage.load_restart_options()
    if options_storage.continue_mode:
        return None, None

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
        options_storage.dataset_yaml_filename = os.path.join(options_storage.output_dir, "input_dataset.yaml")
        pyyaml.dump(dataset_data, open(options_storage.dataset_yaml_filename, 'w'))

    support.check_dataset_reads(dataset_data, options_storage.only_assembler, log)
    if support.dataset_has_only_mate_pairs_libraries(dataset_data):
        support.error('you should specify at least one paired-end or unpaired library (only mate-pairs libraries were found)!')

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

    # dataset section
    cfg["dataset"].__dict__["single_cell"] = options_storage.single_cell
    cfg["dataset"].__dict__["iontorrent"] = options_storage.iontorrent
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

    # assembly
    if not options_storage.only_error_correction:
        if options_storage.k_mers:
            cfg["assembly"].__dict__["iterative_K"] = options_storage.k_mers
        else:
            cfg["assembly"].__dict__["iterative_K"] = options_storage.K_MERS_SHORT
        cfg["assembly"].__dict__["careful"] = options_storage.careful
        cfg["assembly"].__dict__["disable_rr"] = options_storage.disable_rr
        cfg["assembly"].__dict__["diploid_mode"] = options_storage.diploid_mode
        if options_storage.spades_heap_check:
            cfg["assembly"].__dict__["heap_check"] = options_storage.spades_heap_check
        if options_storage.read_buffer_size:
            cfg["assembly"].__dict__["read_buffer_size"] = options_storage.read_buffer_size

    #corrector can work only if contigs exist (not only error correction)
    if (not options_storage.only_error_correction) and options_storage.mismatch_corrector:
        cfg["mismatch_corrector"] = empty_config()
        cfg["mismatch_corrector"].__dict__["skip-masked"] = None
        cfg["mismatch_corrector"].__dict__["bwa"] = os.path.join(bin_home, "bwa-spades")
        cfg["mismatch_corrector"].__dict__["threads"] = options_storage.threads
        cfg["mismatch_corrector"].__dict__["output-dir"] = options_storage.output_dir

    return cfg, dataset_data


def check_cfg_for_restart_from(cfg):
    if options_storage.restart_from == 'ec' and ("error_correction" not in cfg):
        support.error("failed to restart from read error correction because this stage was not specified!")
    if options_storage.restart_from == 'mc' and ("mismatch_corrector" not in cfg):
        support.error("failed to restart from mismatch correction because this stage was not specified!")
    if options_storage.restart_from == 'as' or options_storage.restart_from.startswith('k'):
        if "assembly" not in cfg:
            support.error("failed to restart from assembling because this stage was not specified!")
        if options_storage.restart_from.startswith('k'):
            correct_k = False
            k_to_check = options_storage.k_mers
            if not k_to_check:
                if options_storage.auto_K_allowed():
                    k_to_check = list(set(options_storage.K_MERS_SHORT + options_storage.K_MERS_150 + options_storage.K_MERS_250))
                else:
                    k_to_check = options_storage.K_MERS_SHORT
            for k in k_to_check:
                if options_storage.restart_from == ("k%d" % k) or options_storage.restart_from.startswith("k%d:" % k):
                    correct_k = True
                    break
            if not correct_k:
                k_str = options_storage.restart_from[1:]
                if k_str.find(":") != -1:
                    k_str = k_str[:k_str.find(":")]
                support.error("failed to restart from K=%s because this K was not specified!" % k_str)


def get_options_from_params(params_filename):
    if not os.path.isfile(params_filename):
        return None, None
    params = open(params_filename, 'r')
    cmd_line = params.readline()
    params.close()
    spades_py_pos = cmd_line.find('spades.py')
    if spades_py_pos == -1:
        return None, None
    return cmd_line, cmd_line[spades_py_pos + len('spades.py'):].split()


def main(args):
    os.environ["LC_ALL"] = "C"

    if len(args) == 1:
        options_storage.usage(spades_version)
        sys.exit(0)

    log = logging.getLogger('spades')
    log.setLevel(logging.DEBUG)

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    support.check_binaries(bin_home, log)

    # parse options and safe all parameters to cfg
    cfg, dataset_data = fill_cfg(args, log)

    if options_storage.continue_mode:
        cmd_line, options = get_options_from_params(os.path.join(options_storage.output_dir, "params.txt"))
        if not options:
            support.error("failed to parse command line of the previous run! Please restart from the beginning or specify another output directory.")
        cfg, dataset_data = fill_cfg(options, log)
        if options_storage.restart_from:
            check_cfg_for_restart_from(cfg)
        options_storage.continue_mode = True

    log_filename = os.path.join(cfg["common"].output_dir, "spades.log")
    if options_storage.continue_mode:
        log_handler = logging.FileHandler(log_filename, mode='a')
    else:
        log_handler = logging.FileHandler(log_filename, mode='w')
    log.addHandler(log_handler)

    if options_storage.continue_mode:
        log.info("\n======= SPAdes pipeline continued. Log can be found here: " + log_filename + "\n")
        log.info("Restored from " + cmd_line.strip())
        updated_params = "with updated parameters: "
        flag = False
        for v in args[1:]:
            if v == '-o':
                flag = True
                continue
            if flag:
                flag = False
                continue
            updated_params += " " + v
        log.info(updated_params)

        print_used_values(cfg, log)
    else:
        params_filename = os.path.join(cfg["common"].output_dir, "params.txt")
        params_handler = logging.FileHandler(params_filename, mode='w')
        log.addHandler(params_handler)

        command = "Command line:"
        for v in args:
            command += " " + v
        log.info(command)

        print_used_values(cfg, log)
        log.removeHandler(params_handler)

        log.info("\n======= SPAdes pipeline started. Log can be found here: " + log_filename + "\n")

    # splitting interlaced reads and processing Ns in additional contigs if needed
    if support.dataset_has_interlaced_reads(dataset_data) or support.dataset_has_additional_contigs(dataset_data):
        dir_for_split_reads = os.path.join(options_storage.output_dir, 'split_input')
        if support.dataset_has_interlaced_reads(dataset_data):
            if not os.path.isdir(dir_for_split_reads):
                os.makedirs(dir_for_split_reads)
            dataset_data = support.split_interlaced_reads(dataset_data, dir_for_split_reads, log)
        if support.dataset_has_additional_contigs(dataset_data):
            dataset_data = support.process_Ns_in_additional_contigs(dataset_data, dir_for_split_reads, log)
        options_storage.dataset_yaml_filename = os.path.join(options_storage.output_dir, "input_dataset.yaml")
        pyyaml.dump(dataset_data, open(options_storage.dataset_yaml_filename, 'w'))
        cfg["dataset"].yaml_filename = options_storage.dataset_yaml_filename

    try:
        # copying configs before all computations (to prevent its changing at run time)
        tmp_configs_dir = os.path.join(cfg["common"].output_dir, "configs")
        if os.path.isdir(tmp_configs_dir) and not options_storage.continue_mode:
            shutil.rmtree(tmp_configs_dir)
        if not os.path.isdir(tmp_configs_dir):
            shutil.copytree(os.path.join(spades_home, "configs"), tmp_configs_dir)

        corrected_dataset_yaml_filename = ''
        if "error_correction" in cfg:
            STAGE_NAME = "Read error correction"
            bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])
            corrected_dataset_yaml_filename = os.path.join(bh_cfg.output_dir, "corrected.yaml")
            if os.path.isfile(corrected_dataset_yaml_filename) and options_storage.continue_mode \
                and not options_storage.restart_from == "ec":
                log.info("\n===== Skipping %s (already processed). \n" % STAGE_NAME)
            else:
                support.continue_from_here(log)

                if "HEAPCHECK" in os.environ:
                    del os.environ["HEAPCHECK"]
                if "heap_check" in bh_cfg.__dict__:
                    os.environ["HEAPCHECK"] = bh_cfg.heap_check

                if os.path.exists(bh_cfg.output_dir):
                    shutil.rmtree(bh_cfg.output_dir)
                os.makedirs(bh_cfg.output_dir)

                if support.get_lib_ids_by_type(dataset_data, options_storage.LONG_READS_TYPES):
                    not_used_dataset_data = support.get_libs_by_type(dataset_data, options_storage.LONG_READS_TYPES)
                    to_correct_dataset_data = support.rm_libs_by_type(dataset_data, options_storage.LONG_READS_TYPES)
                    to_correct_dataset_yaml_filename = os.path.join(bh_cfg.output_dir, "to_correct.yaml")
                    pyyaml.dump(to_correct_dataset_data, open(to_correct_dataset_yaml_filename, 'w'))
                    bh_cfg.__dict__["dataset_yaml_filename"] = to_correct_dataset_yaml_filename
                else:
                    not_used_dataset_data = None
                    bh_cfg.__dict__["dataset_yaml_filename"] = cfg["dataset"].yaml_filename

                log.info("\n===== %s started. \n" % STAGE_NAME)
                hammer_logic.run_hammer(corrected_dataset_yaml_filename, tmp_configs_dir, bin_home, bh_cfg, not_used_dataset_data,
                    ext_python_modules_home, log)
                log.info("\n===== %s finished. \n" % STAGE_NAME)

        result_contigs_filename = os.path.join(cfg["common"].output_dir, "contigs.fasta")
        result_scaffolds_filename = os.path.join(cfg["common"].output_dir, "scaffolds.fasta")
        misc_dir = os.path.join(cfg["common"].output_dir, "misc")
        ### if mismatch correction is enabled then result contigs are copied to misc directory
        assembled_contigs_filename = os.path.join(misc_dir, "assembled_contigs.fasta")
        assembled_scaffolds_filename = os.path.join(misc_dir, "assembled_scaffolds.fasta")
        if "assembly" in cfg:
            STAGE_NAME = "Assembling"
            spades_cfg = merge_configs(cfg["assembly"], cfg["common"])
            spades_cfg.__dict__["result_contigs"] = result_contigs_filename
            spades_cfg.__dict__["result_scaffolds"] = result_scaffolds_filename

            if options_storage.continue_mode and (os.path.isfile(spades_cfg.result_contigs)
                                                  or ("mismatch_corrector" in cfg and
                                                      os.path.isfile(assembled_contigs_filename)))\
                and not options_storage.restart_from == 'as' \
                and not (options_storage.restart_from and options_storage.restart_from.startswith('k')):

                log.info("\n===== Skipping %s (already processed). \n" % STAGE_NAME)
                # calculating latest_dir for the next stages
                latest_dir = support.get_latest_dir(os.path.join(spades_cfg.output_dir, "K*"))
                if not latest_dir:
                    support.error("failed to continue the previous run! Please restart from previous stages or from the beginning.", log)
            else:
                old_result_files = [result_contigs_filename, result_scaffolds_filename]
                for format in [".fasta", ".fastg"]:
                    for old_result_file in old_result_files:
                        if os.path.isfile(old_result_file[:-6] + format):
                            os.remove(old_result_file[:-6] + format)

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
                    dataset_file.write("single_cell" + '\t' + process_cfg.bool_to_str(cfg["dataset"].single_cell) + '\n')
                    if os.path.isfile(corrected_dataset_yaml_filename):
                        dataset_file.write("reads" + '\t' + process_cfg.process_spaces(corrected_dataset_yaml_filename) + '\n')
                    else:
                        dataset_file.write("reads" + '\t' + process_cfg.process_spaces(cfg["dataset"].yaml_filename) + '\n')
                    if spades_cfg.developer_mode and "reference" in cfg["dataset"].__dict__:
                        dataset_file.write("reference_genome" + '\t')
                        dataset_file.write(process_cfg.process_spaces(cfg["dataset"].reference) + '\n')
                    dataset_file.close()
                spades_cfg.__dict__["dataset"] = dataset_filename

                latest_dir = spades_logic.run_spades(tmp_configs_dir, bin_home, spades_cfg, dataset_data, ext_python_modules_home, log)

                if os.path.isdir(misc_dir) and not options_storage.continue_mode:
                    shutil.rmtree(misc_dir)
                if not os.path.isdir(misc_dir):
                    os.makedirs(misc_dir)

                if options_storage.continue_mode and options_storage.restart_from and options_storage.restart_from.startswith('k'):
                    k_str = options_storage.restart_from[1:]
                    if k_str.find(":") != -1:
                        k_str = k_str[:k_str.find(":")]
                    support.error("failed to continue from K=%s because this K was not processed in the original run!" % k_str, log)
                log.info("\n===== %s finished. \n" % STAGE_NAME)

            #corrector
            if "mismatch_corrector" in cfg and (os.path.isfile(result_contigs_filename) or
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
                    for format in [".fasta", ".fastg"]:
                        if os.path.isfile(old[:-6] + format):
                            shutil.move(old[:-6] + format, new[:-6] + format)

                if options_storage.continue_mode and os.path.isfile(result_contigs_filename) and \
                    (os.path.isfile(result_scaffolds_filename) or not os.path.isfile(assembled_scaffolds_filename)) \
                    and not options_storage.restart_from == 'mc':
                    log.info("\n===== Skipping %s (already processed). \n" % STAGE_NAME)
                else:
                    if options_storage.restart_from == 'mc':
                        support.continue_from_here(log)

                    log.info("\n===== %s started." % STAGE_NAME)
                    # detecting paired-end library with the largest insert size
                    est_params_data = pyyaml.load(open(os.path.join(latest_dir, "final.lib_data"), 'r'))
                    max_IS_library = None
                    for reads_library in est_params_data:
                        if reads_library['type'] == 'paired-end':
                            if not max_IS_library or float(reads_library["insert size mean"]) > float(max_IS_library["insert size mean"]):
                                max_IS_library = reads_library
                    if not max_IS_library:
                        support.error('Mismatch correction cannot be performed without at least one paired-end library!', log)
                    if not max_IS_library["insert size mean"]:
                        support.warning('Failed to estimate insert size for all paired-end libraries. Starting Mismatch correction'
                                        ' based on the first paired-end library and with default insert size.', log)
                    else:
                        cfg["mismatch_corrector"].__dict__["insert-size"] = round(max_IS_library["insert size mean"])
                    yaml_dirname = os.path.dirname(options_storage.dataset_yaml_filename)
                    cfg["mismatch_corrector"].__dict__["1"] = list(map(lambda x: os.path.join(yaml_dirname, x),
                        max_IS_library['left reads']))
                    cfg["mismatch_corrector"].__dict__["2"] = list(map(lambda x: os.path.join(yaml_dirname, x),
                        max_IS_library['right reads']))
                    #TODO: add reads orientation

                    import corrector
                    corrector_cfg = cfg["mismatch_corrector"]
                    args = []
                    for key, values in corrector_cfg.__dict__.items():
                        if key == "output-dir":
                            continue

                        # for processing list of reads
                        if not isinstance(values, list):
                            values = [values]
                        for value in values:
                            if len(key) == 1:
                                args.append('-' + key)
                            else:
                                args.append('--' + key)
                            if value is not None:
                                args.append(value)

                    # processing contigs and scaffolds (or only contigs)
                    for assembly_type, (corrected, assembled) in to_correct.items():
                        if options_storage.continue_mode and os.path.isfile(corrected):
                            log.info("\n== Skipping processing of " + assembly_type + " (already processed)\n")
                            continue

                        support.continue_from_here(log)
                        log.info("\n== Processing of " + assembly_type + "\n")

                        cur_args = args[:]
                        cur_args += ['-c', assembled]
                        tmp_dir_for_corrector = support.get_tmp_dir(prefix="mis_cor_%s_" % assembly_type)
                        cur_args += ['--output-dir', tmp_dir_for_corrector]

                        # correcting
                        corrector.main(cur_args, ext_python_modules_home, log)

                        result_corrected_filename = os.path.join(tmp_dir_for_corrector, "corrected_contigs.fasta")
                        # moving corrected contigs (scaffolds) to SPAdes output dir
                        if os.path.isfile(result_corrected_filename):
                            shutil.move(result_corrected_filename, corrected)

                        if os.path.isdir(tmp_dir_for_corrector):
                            shutil.rmtree(tmp_dir_for_corrector)

                        assembled_fastg = assembled[:-6] + ".fastg"
                        if os.path.isfile(assembled_fastg):
                            support.create_fastg_from_fasta(corrected, assembled_fastg, log)
                    log.info("\n===== %s finished.\n" % STAGE_NAME)

        if not cfg["common"].developer_mode and os.path.isdir(tmp_configs_dir):
            shutil.rmtree(tmp_configs_dir)

        #log.info("")
        if "error_correction" in cfg and os.path.isdir(os.path.dirname(corrected_dataset_yaml_filename)):
            log.info(" * Corrected reads are in " + os.path.dirname(corrected_dataset_yaml_filename) + "/")
        if "assembly" in cfg and os.path.isfile(result_contigs_filename):
            message = " * Assembled contigs are in " + result_contigs_filename
            if os.path.isfile(result_contigs_filename[:-6] + ".fastg"):
                message += " (" + os.path.basename(result_contigs_filename[:-6] + ".fastg") + ")"
            log.info(message)
        if "assembly" in cfg and os.path.isfile(result_scaffolds_filename):
            message = " * Assembled scaffolds are in " + result_scaffolds_filename
            if os.path.isfile(result_scaffolds_filename[:-6] + ".fastg"):
                message += " (" + os.path.basename(result_scaffolds_filename[:-6] + ".fastg") + ")"
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
    except BaseException: # since python 2.5 system-exiting exceptions (e.g. KeyboardInterrupt) are derived from BaseException
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            support.error("exception caught: %s" % exc_type, log)


if __name__ == '__main__':
    main(sys.argv)
