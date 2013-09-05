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
import bh_logic
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
            log.info("Debug mode turned ON")
        else:
            log.info("Debug mode turned OFF")
    log.info("")

    # dataset
    if "dataset" in cfg:
        log.info("Dataset parameters:")

        if cfg["dataset"].single_cell:
            log.info("  Single-cell mode")
        else:
            log.info("  Multi-cell mode (you should set '--sc' flag if input data"\
                     " was obtained with MDA (single-cell) technology")

        log.info("  Reads:")
        dataset_data = pyyaml.load(open(cfg["dataset"].yaml_filename, 'r'))
        dataset_data = support.relative2abs_paths(dataset_data, os.path.dirname(cfg["dataset"].yaml_filename))
        support.pretty_print_reads(dataset_data, log)

    # error correction
    if "error_correction" in cfg:
        log.info("Read error correction parameters:")
        print_value(cfg, "error_correction", "tmp_dir", "Dir for temp files")
        print_value(cfg, "error_correction", "max_iterations", "Iterations")
        print_value(cfg, "error_correction", "qvoffset", "PHRED offset")

        if cfg["error_correction"].gzip_output:
            log.info("  Corrected reads will be compressed (with gzip)")
        else:
            log.info("  Corrected reads will NOT be compressed (with gzip)")

    # assembly
    if "assembly" in cfg:
        log.info("Assembly parameters:")
        print_value(cfg, "assembly", "iterative_K", "k")

    log.info("Other parameters:")
    print_value(cfg, "common", "max_threads", "Threads")
    print_value(cfg, "common", "max_memory", "Memory limit (in Gb)", "  ")
    log.info("")


def check_binaries(binary_dir, log):
    for binary in ["hammer", "spades", "bwa-spades"]:
        binary_path = os.path.join(binary_dir, binary)
        if not os.path.isfile(binary_path):
            support.error("SPAdes binaries not found: " + binary_path +
                          "\nYou can obtain SPAdes binaries in one of two ways:" +
                          "\n1. Download them from http://spades.bioinf.spbau.ru/release" +
                          str(spades_version).strip() + "/SPAdes-" + str(spades_version).strip() + "-Linux.tar.gz" +
                          "\n2. Build source code with ./spades_compile.sh script", log)


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

    # all parameters are stored here
    cfg = dict()
    # dataset is stored here. We are prepared for up to MAX_LIBS_NUMBER paired-end libs and MAX_LIBS_NUMBER mate-pair libs
    dataset_data = [{} for i in range(options_storage.MAX_LIBS_NUMBER * 2)]

    options_storage.continue_mode = False
    for opt, arg in options:
        if opt == '-o':
            options_storage.output_dir = arg
        elif opt == "--tmp-dir":
            options_storage.tmp_dir = arg
        elif opt == "--reference":
            options_storage.reference = support.check_file_existence(arg, 'reference', log)
        elif opt == "--dataset":
            options_storage.dataset_yaml_filename = support.check_file_existence(arg, 'dataset', log)

        elif opt in options_storage.reads_options:
            support.add_to_dataset(opt, arg, dataset_data)

        elif opt == '-k':
            options_storage.k_mers = list(map(int, arg.split(",")))
            for k in options_storage.k_mers:
                if k > 127:
                    support.error('wrong k value ' + str(k) + ': all k values should be less than 128', log)
                if k % 2 == 0:
                    support.error('wrong k value ' + str(k) + ': all k values should be odd', log)

        elif opt == "--sc":
            options_storage.single_cell = True
        elif opt == "--disable-gzip-output":
            options_storage.disable_gzip_output = True

        elif opt == "--only-error-correction":
            if options_storage.only_assembler:
                support.error('you cannot specify --only-error-correction and --only-assembler simultaneously')
            options_storage.only_error_correction = True
        elif opt == "--only-assembler":
            if options_storage.only_error_correction:
                support.error('you cannot specify --only-error-correction and --only-assembler simultaneously')
            options_storage.only_assembler = True

        elif opt == "--bh-heap-check":
            options_storage.bh_heap_check = arg
        elif opt == "--spades-heap-check":
            options_storage.spades_heap_check = arg

        elif opt == "--continue":
            options_storage.continue_mode = True

        elif opt == '-t' or opt == "--threads":
            options_storage.threads = int(arg)
        elif opt == '-m' or opt == "--memory":
            options_storage.memory = int(arg)
        elif opt == "--phred-offset":
            if int(arg) in [33, 64]:
                options_storage.qvoffset = int(arg)
            else:
                support.error('wrong PHRED quality offset value ' + str(arg) + ': should be either 33 or 64', log)
        elif opt == '-i' or opt == "--iterations":
            options_storage.iterations = int(arg)

        elif opt == "--debug":
            options_storage.developer_mode = True

        elif opt == "--rectangles":
            options_storage.rectangles = True

        #corrector
        elif opt == "--mismatch-correction":
            options_storage.mismatch_corrector = True

        elif opt == "--careful":
            options_storage.mismatch_corrector = True
            options_storage.careful = True

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
        else:
            raise ValueError


    if not options_storage.output_dir:
        support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).", log)
    if not os.path.isdir(options_storage.output_dir):
        if options_storage.continue_mode:
            support.error("the output_dir should exist for --continue!", log)
        os.makedirs(options_storage.output_dir)
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
    if options_storage.rectangles and (len(dataset_data) > 1):
        support.error('rectangle graph algorithm for repeat resolution cannot work with multiple libraries!')

    ### FILLING cfg
    cfg["common"] = empty_config()
    cfg["dataset"] = empty_config()
    if not options_storage.only_assembler:
        cfg["error_correction"] = empty_config()
    if not options_storage.only_error_correction:
        cfg["assembly"] = empty_config()

    # common
    cfg["common"].__dict__["output_dir"] = os.path.abspath(options_storage.output_dir)
    cfg["common"].__dict__["max_threads"] = options_storage.threads
    cfg["common"].__dict__["max_memory"] = options_storage.memory
    cfg["common"].__dict__["developer_mode"] = options_storage.developer_mode

    # dataset section
    cfg["dataset"].__dict__["single_cell"] = options_storage.single_cell
    cfg["dataset"].__dict__["yaml_filename"] = os.path.abspath(options_storage.dataset_yaml_filename)
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
        if options_storage.tmp_dir:
            cfg["error_correction"].__dict__["tmp_dir"] = options_storage.tmp_dir
        else:
            cfg["error_correction"].__dict__["tmp_dir"] = cfg["error_correction"].output_dir
        cfg["error_correction"].tmp_dir = os.path.join(os.path.abspath(cfg["error_correction"].tmp_dir), 'tmp')

    # assembly
    if not options_storage.only_error_correction:
        if options_storage.k_mers:
            cfg["assembly"].__dict__["iterative_K"] = options_storage.k_mers
        else:
            cfg["assembly"].__dict__["iterative_K"] = options_storage.k_mers_short
        cfg["assembly"].__dict__["careful"] = options_storage.careful
        if options_storage.spades_heap_check:
            cfg["assembly"].__dict__["heap_check"] = options_storage.spades_heap_check

    #corrector can work only if contigs exist (not only error correction)
    if (not options_storage.only_error_correction) and options_storage.mismatch_corrector:
        cfg["mismatch_corrector"] = empty_config()
        cfg["mismatch_corrector"].__dict__["skip-masked"] = ""
        cfg["mismatch_corrector"].__dict__["bwa"] = os.path.join(bin_home, "bwa-spades")
        cfg["mismatch_corrector"].__dict__["threads"] = options_storage.threads
        cfg["mismatch_corrector"].__dict__["output-dir"] = options_storage.output_dir

    return cfg, dataset_data


def get_options_from_params(params_filename):
    if not os.path.isfile(params_filename):
        return None
    params = open(params_filename, 'r')
    cmd_line = params.readline()
    params.close()
    spades_py_pos = cmd_line.find('spades.py')
    if spades_py_pos == -1:
        return None
    return cmd_line, cmd_line[spades_py_pos + len('spades.py'):].split()


def main():
    os.environ["LC_ALL"] = "C"

    if len(sys.argv) == 1:
        options_storage.usage(spades_version)
        sys.exit(0)

    log = logging.getLogger('spades')
    log.setLevel(logging.DEBUG)

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    check_binaries(bin_home, log)

    # parse options and safe all parameters to cfg
    cfg, dataset_data = fill_cfg(sys.argv, log)

    if options_storage.continue_mode:
        cmd_line, options = get_options_from_params(os.path.join(options_storage.output_dir, "params.txt"))
        if not options:
            support.error("failed to parse command line of the previous run! Please restart from the beginning.")
        cfg, dataset_data = fill_cfg(options, log)
        options_storage.continue_mode = True

    log_filename = os.path.join(cfg["common"].output_dir, "spades.log")
    if options_storage.continue_mode:
        log_handler = logging.FileHandler(log_filename, mode='a')
    else:
        log_handler = logging.FileHandler(log_filename, mode='w')
    log.addHandler(log_handler)

    if options_storage.continue_mode:
        log.info("\n======= SPAdes pipeline continued. Log can be found here: " + log_filename + "\n")
        log.info("Restored from " + cmd_line)
    else:
        params_filename = os.path.join(cfg["common"].output_dir, "params.txt")
        params_handler = logging.FileHandler(params_filename, mode='w')
        log.addHandler(params_handler)

        command = "Command line:"
        for v in sys.argv:
            command += " " + v
        log.info(command)

        print_used_values(cfg, log)
        log.removeHandler(params_handler)

        log.info("\n======= SPAdes pipeline started. Log can be found here: " + log_filename + "\n")

    # splitting interlaced reads if needed
    if support.dataset_has_interlaced_reads(dataset_data):
        dir_for_split_reads = os.path.join(os.path.abspath(options_storage.output_dir), 'split_reads')
        if not os.path.isdir(dir_for_split_reads):
            os.makedirs(dir_for_split_reads)
        dataset_data = support.split_interlaced_reads(dataset_data, dir_for_split_reads, log)
        options_storage.dataset_yaml_filename = os.path.join(options_storage.output_dir, "input_dataset.yaml")
        pyyaml.dump(dataset_data, open(options_storage.dataset_yaml_filename, 'w'))
        cfg["dataset"].yaml_filename = os.path.abspath(options_storage.dataset_yaml_filename)

    try:
        # copying configs before all computations (to prevent its changing at run time)
        tmp_configs_dir = os.path.join(cfg["common"].output_dir, "configs")
        if os.path.isdir(tmp_configs_dir) and not options_storage.continue_mode:
            shutil.rmtree(tmp_configs_dir)
        if not os.path.isdir(tmp_configs_dir):
            shutil.copytree(os.path.join(spades_home, "configs"), tmp_configs_dir)

        corrected_dataset_yaml_filename = ''
        if "error_correction" in cfg:
            bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])
            bh_cfg.__dict__["dataset_yaml_filename"] = cfg["dataset"].yaml_filename
            corrected_dataset_yaml_filename = os.path.join(bh_cfg.output_dir, "corrected.yaml")
            if os.path.isfile(corrected_dataset_yaml_filename) and options_storage.continue_mode:
                log.info("\n===== Skipping read error correction (already processed). \n")
            else:
                options_storage.continue_mode = False # continue from here

                if "HEAPCHECK" in os.environ:
                    del os.environ["HEAPCHECK"]
                if "heap_check" in bh_cfg.__dict__:
                    os.environ["HEAPCHECK"] = bh_cfg.heap_check

                if os.path.exists(bh_cfg.output_dir):
                    shutil.rmtree(bh_cfg.output_dir)

                os.makedirs(bh_cfg.output_dir)
                if not os.path.exists(bh_cfg.tmp_dir):
                    os.makedirs(bh_cfg.tmp_dir)

                log.info("\n===== Read error correction started. \n")
                bh_logic.run_bh(corrected_dataset_yaml_filename, tmp_configs_dir, bin_home, bh_cfg,
                    ext_python_modules_home, log)
                log.info("\n===== Read error correction finished. \n")

        result_contigs_filename = os.path.join(cfg["common"].output_dir, "contigs.fasta")
        result_scaffolds_filename = os.path.join(cfg["common"].output_dir, "scaffolds.fasta")
        misc_dir = os.path.join(cfg["common"].output_dir, "misc")
        ### if mismatch correction is enabled then result contigs are copied to misc directory
        assembled_contigs_filename = os.path.join(misc_dir, "assembled_contigs.fasta")
        assembled_scaffolds_filename = os.path.join(misc_dir, "assembled_scaffolds.fasta")
        if "assembly" in cfg:
            spades_cfg = merge_configs(cfg["assembly"], cfg["common"])
            spades_cfg.__dict__["result_contigs"] = result_contigs_filename
            spades_cfg.__dict__["result_scaffolds"] = result_scaffolds_filename
            spades_cfg.__dict__["additional_contigs"] = os.path.join(spades_cfg.output_dir, "simplified_contigs.fasta")

            if options_storage.continue_mode and (os.path.isfile(spades_cfg.result_contigs)
                                                  or ("mismatch_corrector" in cfg and
                                                      os.path.isfile(assembled_contigs_filename))):
                log.info("\n===== Skipping assembling (already processed). \n")
                # calculating latest_dir for the next stages
                latest_dir = support.get_latest_dir(os.path.join(spades_cfg.output_dir, "K*"))
                if not latest_dir:
                    support.error("failed to continue the previous run! Please restart from the beginning.")
            else:
                if os.path.isfile(corrected_dataset_yaml_filename):
                    dataset_data = pyyaml.load(open(corrected_dataset_yaml_filename, 'r'))
                    dataset_data = support.relative2abs_paths(dataset_data, os.path.dirname(corrected_dataset_yaml_filename))
                if support.dataset_has_paired_reads(dataset_data):
                    spades_cfg.__dict__["paired_mode"] = True
                else:
                    spades_cfg.__dict__["paired_mode"] = False

                if options_storage.rectangles:
                    spades_cfg.__dict__["resolving_mode"] = "rectangles"

                if "HEAPCHECK" in os.environ:
                    del os.environ["HEAPCHECK"]
                if "heap_check" in spades_cfg.__dict__:
                    os.environ["HEAPCHECK"] = spades_cfg.heap_check

                log.info("\n===== Assembling started.\n")

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
                        dataset_file.write(process_cfg.process_spaces(os.path.abspath(cfg["dataset"].reference)) + '\n')
                    dataset_file.close()
                spades_cfg.__dict__["dataset"] = dataset_filename

                latest_dir = spades_logic.run_spades(tmp_configs_dir, bin_home, spades_cfg, log)

                #rectangles
                if spades_cfg.paired_mode and options_storage.rectangles:
                    if options_storage.continue_mode: # TODO: continue mode
                        support.warning("sorry, --continue doesn't work with --rectangles yet. Skipping repeat resolving.")
                    else:
                        sys.path.append(os.path.join(python_modules_home, "rectangles"))
                        import rrr

                        rrr_input_dir = os.path.join(latest_dir, "saves")
                        rrr_outpath = os.path.join(spades_cfg.output_dir, "rectangles")
                        if not os.path.exists(rrr_outpath):
                            os.mkdir(rrr_outpath)

                        rrr_reference_information_file = os.path.join(rrr_input_dir,
                            "late_pair_info_counted_etalon_distance.txt")
                        rrr_test_util = rrr.TestUtils(rrr_reference_information_file,
                            os.path.join(rrr_outpath, "rectangles.log"))
                        rrr.resolve(rrr_input_dir, rrr_outpath, rrr_test_util, "", cfg["dataset"].single_cell, spades_cfg.careful)

                        shutil.copyfile(os.path.join(rrr_outpath, "rectangles_extend_before_scaffold.fasta"), spades_cfg.result_contigs)
                        shutil.copyfile(os.path.join(rrr_outpath, "rectangles_extend.fasta"), spades_cfg.result_scaffolds)

                        if not spades_cfg.developer_mode:
                            if os.path.exists(rrr_input_dir):
                                shutil.rmtree(rrr_input_dir)
                            if os.path.exists(rrr_outpath):
                                shutil.rmtree(rrr_outpath, True)
                            if os.path.exists(rrr_outpath):
                                os.system('rm -r ' + rrr_outpath)
                                #EOR

                if os.path.isdir(misc_dir) and not options_storage.continue_mode:
                    shutil.rmtree(misc_dir)
                if not os.path.isdir(misc_dir):
                    os.makedirs(misc_dir)
                    if os.path.isfile(spades_cfg.additional_contigs):
                        shutil.move(spades_cfg.additional_contigs, misc_dir)

                log.info("\n===== Assembling finished. \n")

            #corrector
            if "mismatch_corrector" in cfg and (os.path.isfile(result_contigs_filename) or
                                                (options_storage.continue_mode and os.path.isfile(assembled_contigs_filename))):
                to_correct = dict()
                to_correct["contigs"] = (result_contigs_filename, assembled_contigs_filename)
                if os.path.isfile(result_scaffolds_filename) or (options_storage.continue_mode and
                                                                 os.path.isfile(assembled_scaffolds_filename)):
                    to_correct["scaffolds"] = (result_scaffolds_filename, assembled_scaffolds_filename)

                # moving assembled contigs (scaffolds) to misc dir
                for k, (old, new) in to_correct.items():
                    if options_storage.continue_mode and os.path.isfile(new):
                        continue
                    shutil.move(old, new)

                if options_storage.continue_mode and os.path.isfile(result_contigs_filename) and \
                    (os.path.isfile(result_scaffolds_filename) or not os.path.isfile(assembled_scaffolds_filename)):
                    log.info("\n===== Skipping mismatch correction (already processed). \n")
                else:
                    log.info("\n===== Mismatch correction started.")

                    # detecting paired-end library with the largest insert size
                    dataset_data = pyyaml.load(open(options_storage.dataset_yaml_filename, 'r')) ### initial dataset, i.e. before error correction
                    dataset_data = support.relative2abs_paths(dataset_data, os.path.dirname(options_storage.dataset_yaml_filename))
                    paired_end_libraries_ids = []
                    for id, reads_library in enumerate(dataset_data):
                        if reads_library['type'] == 'paired-end':
                            paired_end_libraries_ids.append(id)
                    if not len(paired_end_libraries_ids):
                        support.error('Mismatch correction cannot be performed without at least one paired-end library!')
                    estimated_params = load_config_from_file(os.path.join(latest_dir, "_est_params.info"))
                    max_insert_size = -1
                    target_paired_end_library_id = -1
                    for id in paired_end_libraries_ids:
                        if float(estimated_params.__dict__["insert_size_" + str(id)]) > max_insert_size:
                            max_insert_size = float(estimated_params.__dict__["insert_size_" + str(id)])
                            target_paired_end_library_id = id
                    yaml_dirname = os.path.dirname(options_storage.dataset_yaml_filename)
                    cfg["mismatch_corrector"].__dict__["1"] = list(map(lambda x: os.path.join(yaml_dirname, x),
                        dataset_data[target_paired_end_library_id]['left reads']))
                    cfg["mismatch_corrector"].__dict__["2"] = list(map(lambda x: os.path.join(yaml_dirname, x),
                        dataset_data[target_paired_end_library_id]['right reads']))
                    cfg["mismatch_corrector"].__dict__["insert-size"] = round(max_insert_size)
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
                            if value:
                                args.append(value)

                    # processing contigs and scaffolds (or only contigs)
                    for k, (corrected, assembled) in to_correct.items():
                        if options_storage.continue_mode and os.path.isfile(corrected):
                            log.info("\n== Skipping processing of " + k + " (already processed)\n")
                            continue

                        options_storage.continue_mode = False
                        log.info("\n== Processing of " + k + "\n")

                        cur_args = args[:]
                        cur_args += ['-c', assembled]
                        tmp_dir_for_corrector = os.path.join(corrector_cfg.__dict__["output-dir"], "mismatch_corrector_" + k)
                        cur_args += ['--output-dir', tmp_dir_for_corrector]

                        # correcting
                        corrector.main(cur_args, ext_python_modules_home, log)

                        result_corrected_filename = os.path.abspath(os.path.join(tmp_dir_for_corrector, "corrected_contigs.fasta"))
                        # moving corrected contigs (scaffolds) to SPAdes output dir
                        if os.path.isfile(result_corrected_filename):
                            shutil.move(result_corrected_filename, corrected)

                        if os.path.isdir(tmp_dir_for_corrector):
                            shutil.rmtree(tmp_dir_for_corrector)

                    log.info("\n===== Mismatch correction finished.\n")

        if not cfg["common"].developer_mode and os.path.isdir(tmp_configs_dir):
            shutil.rmtree(tmp_configs_dir)

        #log.info("")
        if os.path.isdir(os.path.dirname(corrected_dataset_yaml_filename)):
            log.info(" * Corrected reads are in " + os.path.dirname(corrected_dataset_yaml_filename) + "/")
        if os.path.isfile(result_contigs_filename):
            log.info(" * Assembled contigs are in " + result_contigs_filename)
        if os.path.isfile(result_scaffolds_filename):
            log.info(" * Assembled scaffolds are in " + result_scaffolds_filename)
        #log.info("")

        #breaking scaffolds
        if os.path.isfile(result_scaffolds_filename):
            if not os.path.isdir(misc_dir):
                os.makedirs(misc_dir)
            result_broken_scaffolds = os.path.join(misc_dir, "broken_scaffolds.fasta")
            threshold = 3
            if not os.path.isfile(result_broken_scaffolds) or not options_storage.continue_mode:
                support.break_scaffolds(result_scaffolds_filename, threshold, result_broken_scaffolds)
                #log.info(" * Scaffolds broken by " + str(threshold) + " Ns are in " + result_broken_scaffolds)

        ### printing WARNINGS SUMMARY
        if not support.log_warnings(log):
            log.info("\n======= SPAdes pipeline finished.")  # otherwise it finished WITH WARNINGS

        log.info("\nSPAdes log can be found here: " + log_filename)
        log.info("")
        log.info("Thank you for using SPAdes!")
        log.removeHandler(log_handler)

    except Exception:
        _, exc, _ = sys.exc_info()
        log.exception(exc)
        support.error("exception caught", log)


if __name__ == '__main__':
    main()
