#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import logging
import platform

import spades_init

spades_init.init()
spades_home = spades_init.spades_home
execution_home = os.path.join(spades_home, 'bin')

import support
from process_cfg import *
import bh_logic
import spades_logic


def print_used_values(cfg, log):
    def print_value(cfg, section, param, pretty_param="", margin="  "):
        if not pretty_param:
            pretty_param = param.capitalize().replace('_', ' ')
        line = margin + pretty_param
        if param in cfg[section].__dict__:
            line +=  ": " + str(cfg[section].__dict__[param])
        else:
            if param.find("offset") != -1:
                line += " will be auto-detected"
        log.info(line)

    log.info("")

    # system info
    log.info("System information:")
    try:
        log.info("  OS: " + platform.platform())
        # for more deatils: '[' + str(platform.uname()) + ']'
        log.info("  Python version: " + str(sys.version_info[0]) + "." + str(sys.version_info[1]) + '.' + str(sys.version_info[2]))
        # for more details: '[' + str(sys.version_info) + ']'
    except:
        log.info("  Problem occured when getting system information")
    log.info("")

    # main
    print_value(cfg, "common", "output_dir", "", "")
    if ("error_correction" in cfg) and (not "assembly" in cfg):
        log.info("Mode: ONLY error correction (without assembler)")
    elif (not "error_correction" in cfg) and ("assembly" in cfg):
        log.info("Mode: ONLY assembler (without error correction)")
    else:
        log.info("Mode: error correction and assembler")
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

        no_single = True
        no_paired = True
        for k, v in cfg["dataset"].__dict__.items():
            if k.startswith("single_reads") or k.startswith("paired_reads"):
                if k.startswith("paired_reads"):
                    no_paired = False
                    if not isinstance(v, list):
                        log.info("  Paired reads (file with interlaced reads):")
                    else:
                        log.info("  Paired reads (files with left and right reads):")
                else:
                    no_single = False
                    log.info("  Single reads:")
                if not isinstance(v, list):
                    v = [v]
                for reads_file in v:
                    log.info("    "  + os.path.abspath(os.path.expandvars(reads_file)))

        if no_paired:
            log.info("  Paired reads was not specified")
        if no_single:
            log.info("  Single reads was not specified")

    # error correction
    if "error_correction" in cfg:
        log.info("Error correction parameters:")
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

        if cfg["assembly"].generate_sam_files:
            log.info("  SAM file will be generated")
        else:
            log.info("  SAM file will NOT be generated (WARNING: SAM file are required "\
                   "for some of postprocessing tools)")

        if cfg["assembly"].gap_closer:
            log.info("  The gap closer will be used")
        else:
            log.info("  The gap closer will NOT be used")


    log.info("Other parameters:")
    print_value(cfg, "common", "max_threads", "Threads")
    print_value(cfg, "common", "max_memory", "Memory limit (in Gb)", "  ")
    log.info("")


def check_config(cfg):
    ## special case: 0 iterations for the error correction means "No error correction!"

    if ("error_correction" in cfg) and ("max_iterations" in cfg["error_correction"].__dict__):
        if cfg["error_correction"].max_iterations == 0:
            del cfg["error_correction"]

    ## checking mandatory sections

    if (not "dataset" in cfg):
        support.error("wrong config! You should specify 'dataset' section!")
        return False

    if (not "error_correction" in cfg) and (not "assembly" in cfg):
        support.error("wrong options! You should specify either '--only-error-correction' (for reads"
              " error correction) or '--only-assembler' (for assembling) or none of these options (for both)!")
        return False

    ## checking existence of all files in dataset section

    no_files_with_reads = True
    for k, v in cfg["dataset"].__dict__.items():
        if k.startswith("single_reads") or k.startswith("paired_reads"):
            no_files_with_reads = False
            if not isinstance(v, list):
                v = [v]
            for reads_file in v:
                if not os.path.isfile(os.path.expandvars(reads_file)):
                    support.error("file with reads doesn't exist! " + os.path.expandvars(reads_file))
                    return False
                else:
                    ext = os.path.splitext(os.path.expandvars(reads_file))[1]
                    if ext not in ['.fa', '.fasta', '.fq', '.fastq', '.gz']:
                        support.error("file with reads has unsupported format (only .fa, .fasta, .fq,"
                              " .fastq, .gz are supported)! " + os.path.expandvars(reads_file))
                        return False

    if no_files_with_reads:
        support.error("wrong options! You should specify at least one file with reads!")
        return False

    ## setting default values if needed

    # common 
    if not "output_dir" in cfg["common"].__dict__:
        cfg["common"].__dict__["output_dir"] = 'spades_output'

    if not "max_threads" in cfg["common"].__dict__:
        cfg["common"].__dict__["max_threads"] = 16

    if not "max_memory" in cfg["common"].__dict__:
        cfg["common"].__dict__["max_memory"] = 250

    if not "output_to_console" in cfg["common"].__dict__:
        cfg["common"].__dict__["output_to_console"] = True

    if not "developer_mode" in cfg["common"].__dict__:
        cfg["common"].__dict__["developer_mode"] = False

    cfg["common"].output_dir = os.path.abspath(os.path.expandvars(cfg["common"].output_dir))

    # dataset
    if "dataset" in cfg:
        if not "single_cell" in cfg["dataset"].__dict__:
            cfg["dataset"].__dict__["single_cell"] = False        

    # error_correction
    if "error_correction" in cfg:
        if not "max_iterations" in cfg["error_correction"].__dict__:
            cfg["error_correction"].__dict__["max_iterations"] = 1
        if not "gzip_output" in cfg["error_correction"].__dict__:
            cfg["error_correction"].__dict__["gzip_output"] = True
        if not "tmp_dir" in cfg["error_correction"].__dict__:
            cfg["error_correction"].__dict__["tmp_dir"] = os.path.join(
                cfg["common"].output_dir, os.path.join('corrected', 'tmp'))
        cfg["error_correction"].tmp_dir = os.path.abspath(
            os.path.expandvars(cfg["error_correction"].tmp_dir))

    # assembly
    if "assembly" in cfg:
        if not "iterative_K" in cfg["assembly"].__dict__:
            cfg["assembly"].__dict__["iterative_K"] = [21, 33, 55]
        if not "generate_sam_files" in cfg["assembly"].__dict__:
            cfg["assembly"].__dict__["generate_sam_files"] = False
        if not "gap_closer" in cfg["assembly"].__dict__:
            cfg["assembly"].__dict__["gap_closer"] = True

    return True


def check_binaries(binary_dir):
    for binary in ["hammer", "spades"]:
        binary_path = os.path.join(binary_dir, binary)
        if not os.path.isfile(binary_path):
            support.error("SPAdes binary file not found: " + binary_path +
                  "\nYou can obtain SPAdes binaries in one of two ways:" + 
                  "\n1. Download the binaries from SPAdes server with ./spades_download_binary.py script" + 
                  "\n2. Build source code with ./spades_compile.py script")
            return False
    return True


# for spades.py --test
def fill_test_config(cfg):
    # common
    cfg["common"] = load_config_from_vars(dict())
    cfg["common"].__dict__["output_dir"] = 'spades_test'

    # dataset
    cfg["dataset"] = load_config_from_vars(dict())
    cfg["dataset"].__dict__["paired_reads"] = \
        [os.path.join(spades_home, "test_dataset/ecoli_1K_1.fq.gz"), 
         os.path.join(spades_home, "test_dataset/ecoli_1K_2.fq.gz")]
    cfg["dataset"].__dict__["single_cell"] = False

    # error_correction (load default params)
    cfg["error_correction"] = load_config_from_vars(dict())    

    # assembly (load default params)
    cfg["assembly"] = load_config_from_vars(dict())


long_options = "12= threads= memory= tmp-dir= iterations= phred-offset= sc "\
               "generate-sam-file only-error-correction only-assembler "\
               "disable-gap-closer disable-gzip-output help test debug reference= "\
               "bh-heap-check= spades-heap-check= help-hidden "\
               "config-file= use-jemalloc dataset=".split()
short_options = "o:1:2:s:k:t:m:i:h"


def check_file(f, message=''):
    if not os.path.isfile(f):
        support.error("file not found (%s): %s" % (message, f))
    return f


def usage(show_hidden=False):
    print >> sys.stderr, "SPAdes genome assembler"
    print >> sys.stderr, "Usage:", sys.argv[0], "[options] -o <output_dir>"
    print >> sys.stderr, ""
    print >> sys.stderr, "Options:"
    print >> sys.stderr, "-o\t<output_dir>\tdirectory to store all the resulting files (required)"
    print >> sys.stderr, "--sc\t\t\tthis flag is required for MDA (single-cell)"\
                         " data"
    print >> sys.stderr, "--12\t<filename>\tfile with interlaced left and right"\
                         " paired end reads"
    print >> sys.stderr, "-1\t<filename>\tfile with left paired end reads"
    print >> sys.stderr, "-2\t<filename>\tfile with right paired end reads"
    print >> sys.stderr, "-s\t<filename>\tfile with unpaired reads"
    print >> sys.stderr, "--generate-sam-file\tforces SPAdes to generate SAM-file"

    print >> sys.stderr, ""
    print >> sys.stderr, "Advanced options:"
    print >> sys.stderr, "-t/--threads\t<int>\t\tnumber of threads [default: 16]"
    print >> sys.stderr, "-m/--memory\t<int>\t\tRAM limit for SPAdes in Gb"\
                         " (terminates if exceeded) [default: 250]"
    print >> sys.stderr, "--tmp-dir\t<dirname>\tdirectory for error correction's"\
                         " temp files"
    print >> sys.stderr, "\t\t\t\t[default: <output_dir>/corrected/tmp]"
    print >> sys.stderr, "-k\t<int,int,...>\t\tcomma-separated list of k-mer sizes"\
                         " (must be odd and less than 128)"
    print >> sys.stderr, "\t\t\t\t[default: 21,33,55]"
    print >> sys.stderr, "-i/--iterations\t<int>\t\tnumber of iterations for error"\
                         " correction"
    print >> sys.stderr, "--phred-offset\t<33 or 64>\tPHRED quality offset in the"\
                         " input reads (33 or 64) [default: auto-detect]"
    print >> sys.stderr, "--only-error-correction\t\trun only error correction"\
                         " (without assembler)"
    print >> sys.stderr, "--only-assembler\t\trun only assembler (without error"\
                         " correction)"
    print >> sys.stderr, "--disable-gap-closer\t\tforces SPAdes not to use the gap"\
                         " closer"
    print >> sys.stderr, "--disable-gzip-output\t\tforces error correction not to"\
                         " compress the corrected reads" 

    if show_hidden:
        print >> sys.stderr, ""
        print >> sys.stderr, "HIDDEN options:"  
        print >> sys.stderr, "--config-file\t<filename>\tconfiguration file for spades.py"\
                             " (WARN: all other options will be skipped)"   
        print >> sys.stderr, "--dataset\t<filename>\tfile with dataset description"\
                             " (WARN: works exclusively in --only-assembler mode)" 
        print >> sys.stderr, "--reference\t<filename>\tfile with reference for deep analysis"\
                             " (only in debug mode)" 
        print >> sys.stderr, "--bh-heap-check\t\t<value>\tset HEAPCHECK environment variable"\
                             " for BayesHammer"
        print >> sys.stderr, "--spades-heap-check\t<value>\tset HEAPCHECK environment variable"\
                             " for SPAdes"
        print >> sys.stderr, "--use-jemalloc\t\tall binaries are run with jemalloc.sh"

    print >> sys.stderr, ""
    print >> sys.stderr, "--test\t\trun SPAdes on toy dataset"
    print >> sys.stderr, "--debug\t\trun SPAdes in debug mode"
    print >> sys.stderr, "-h/--help\tprint this usage message"
    if show_hidden:
        print >> sys.stderr, "--help-hidden\tprint this usage message with all hidden options"


def main():
    os.environ["LC_ALL"] = "C"
   
    CONFIG_FILE = ""
    options = None
    TEST = False

    if len(sys.argv) == 1:
        usage()
        sys.exit(0)
        
    try:
        options, not_options = getopt.gnu_getopt(sys.argv, short_options, long_options)
    except (getopt.GetoptError, err):
        print >> sys.stderr, err
        print >> sys.stderr
        usage()
        sys.exit(1)

    # all parameters are stored here
    cfg = dict()

    if options:
        output_dir = ""
        tmp_dir = ""
        reference = ""
        dataset = ""

        paired = []
        paired1 = []
        paired2 = []
        single = []
        k_mers = []

        single_cell = False
        disable_gzip_output = False
        generate_sam_files = False
        disable_gap_closer = False

        only_error_correction = False
        only_assembler = False

        bh_heap_check = ""
        spades_heap_check = ""

        developer_mode = False

        use_jemalloc = False

        threads = None
        memory = None
        qvoffset = None
        iterations = None

        for opt, arg in options:
            if opt == '-o':
                output_dir = arg
            elif opt == "--tmp-dir":
                tmp_dir = arg
            elif opt == "--reference":
                reference = check_file(arg, 'reference')
            elif opt == "--dataset":
                dataset = check_file(arg, 'dataset')

            elif opt == "--12":
                paired.append(check_file(arg, 'paired'))
            elif opt == '-1':
                paired1.append(check_file(arg, 'left paired'))
            elif opt == '-2':
                paired2.append(check_file(arg, 'right paired'))
            elif opt == '-s':
                single.append(check_file(arg, 'single'))
            elif opt == '-k':
                k_mers = map(int, arg.split(","))
                for k in k_mers:
                    if k > 127:
                        support.error('wrong k value ' + str(k) + ': all k values should be less than 128')
                    if k % 2 == 0:
                        support.error('wrong k value ' + str(k) + ': all k values should be odd')

            elif opt == "--sc":
                single_cell = True
            elif opt == "--disable-gzip-output":
                disable_gzip_output = True
            elif opt == "--generate-sam-file":
                generate_sam_files = True
            elif opt == "--disable-gap-closer":
                disable_gap_closer = True

            elif opt == "--only-error-correction":
                only_error_correction = True
            elif opt == "--only-assembler":
                only_assembler = True

            elif opt == "--bh-heap-check":
                bh_heap_check = arg
            elif opt == "--spades-heap-check":
                spades_heap_check = arg
            
            elif opt == '-t' or opt == "--threads":
                threads = int(arg)
            elif opt == '-m' or opt == "--memory":
                memory = int(arg)
            elif opt == "--phred-offset":
                qvoffset = int(arg)
            elif opt == '-i' or opt == "--iterations":
                iterations = int(arg)

            elif opt == "--debug":
                developer_mode = True

            elif opt == "--use-jemalloc":
                use_jemalloc = True

            elif opt == '-h' or opt == "--help":
                usage()
                sys.exit(0)
            elif opt == "--help-hidden":
                usage(True)
                sys.exit(0)

            elif opt == "--config-file":
                CONFIG_FILE = check_file(arg, 'config file')
                break

            elif opt == "--test": 
                fill_test_config(cfg)
                TEST = True
                break
            else:
                raise ValueError

        if not CONFIG_FILE and not TEST:
            if not output_dir:
                support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).")

            if len(paired1) != len(paired2):
                support.error("the number of files with left paired reads is not equal to the"
                      " number of files with right paired reads!")
    
            # processing hidden option "--dataset"
            if dataset: 
                if not only_assembler:
                    support.error("hidden option --dataset works exclusively in --only-assembler mode!")
                # reading info about dataset from provided dataset file
                cfg["dataset"] = load_config_from_file(dataset)
                # correcting reads relative pathes (reads can be specified relatively to provided dataset file)
                for k, v in cfg["dataset"].__dict__.iteritems():
                    if k.find("_reads") != -1:
                        corrected_reads_filenames = []
                        if not isinstance(v, list):
                            v = [v]
                        for reads_filename in v:
                            corrected_reads_filename = os.path.join(os.path.dirname(dataset), reads_filename)
                            check_file(corrected_reads_filename, k + " in " + dataset)
                            corrected_reads_filenames.append(corrected_reads_filename)
                        cfg["dataset"].__dict__[k] = corrected_reads_filenames
            else:
                if not paired and not paired1 and not single:
                    support.error("you should specify either paired reads (-1, -2 or -12) or single reads (-s) or both!")

                # creating empty "dataset" section
                cfg["dataset"] = load_config_from_vars(dict())

                # filling reads
                paired_counter = 0
                if paired:
                    for read in paired:
                        paired_counter += 1
                        cfg["dataset"].__dict__["paired_reads#" + str(paired_counter)] = read

                if paired1:
                    for i in range(len(paired1)):
                        paired_counter += 1
                        cfg["dataset"].__dict__["paired_reads#" + str(paired_counter)] = [paired1[i],
                                                                                          paired2[i]]
                        #cfg["dataset"].__dict__["paired_reads"] = [paired1[i], paired2[i]]

                if single:
                    cfg["dataset"].__dict__["single_reads"] = single

                # filling other parameters
                cfg["dataset"].__dict__["single_cell"] = single_cell
                if developer_mode and reference:
                    cfg["dataset"].__dict__["reference"] = reference

            # filling cfg
            cfg["common"] = load_config_from_vars(dict())
            if not only_assembler:
                cfg["error_correction"] = load_config_from_vars(dict())
            if not only_error_correction:
                cfg["assembly"] = load_config_from_vars(dict())

            # common
            if output_dir:
                cfg["common"].__dict__["output_dir"] = output_dir
            if threads:
                cfg["common"].__dict__["max_threads"] = threads
            if memory:
                cfg["common"].__dict__["max_memory"] = memory
            if developer_mode:
                cfg["common"].__dict__["developer_mode"] = developer_mode
            if use_jemalloc:
                cfg["common"].__dict__["use_jemalloc"] = use_jemalloc

            # error correction
            if not only_assembler:
                if tmp_dir:
                    cfg["error_correction"].__dict__["tmp_dir"] = tmp_dir
                if qvoffset:
                    cfg["error_correction"].__dict__["qvoffset"] = qvoffset
                if iterations != None:
                    cfg["error_correction"].__dict__["max_iterations"] = iterations
                if bh_heap_check:
                    cfg["error_correction"].__dict__["heap_check"] = bh_heap_check
                cfg["error_correction"].__dict__["gzip_output"] = not disable_gzip_output

            # assembly
            if not only_error_correction:
                if k_mers:
                    cfg["assembly"].__dict__["iterative_K"] = k_mers
                if spades_heap_check:
                    cfg["assembly"].__dict__["heap_check"] = spades_heap_check
                cfg["assembly"].__dict__["generate_sam_files"] = generate_sam_files
                cfg["assembly"].__dict__["gap_closer"] = not disable_gap_closer

    if CONFIG_FILE:
        cfg = load_config_from_info_file(CONFIG_FILE)
        os.environ["cfg"] = os.path.dirname(os.path.abspath(CONFIG_FILE))

    if not CONFIG_FILE and not options and not TEST:
        usage()
        sys.exit(1)

    if not check_config(cfg):
        return
        
    if not check_binaries(execution_home):
        return

    print("\n======= SPAdes pipeline started\n")

    if not os.path.isdir(cfg["common"].output_dir):
        os.makedirs(cfg["common"].output_dir)

    log = logging.getLogger('params')
    log.setLevel(logging.DEBUG)

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    params_filename = os.path.join(cfg["common"].output_dir, "params.txt")
    params = logging.FileHandler(params_filename, mode='w')
    log.addHandler(params)

    if CONFIG_FILE:
        log.info("Using config file: " + CONFIG_FILE)
    else:
        command = "Command line:"
        for v in sys.argv:
            command += " " + v
        log.info(command)

    print_used_values(cfg, log)

    bh_dataset_filename = ""
    if "error_correction" in cfg:
        bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])

        if "HEAPCHECK" in os.environ:
            del os.environ["HEAPCHECK"]
        if "heap_check" in bh_cfg.__dict__:
            os.environ["HEAPCHECK"] = bh_cfg.heap_check

        bh_cfg.output_dir = os.path.join(os.path.expandvars(bh_cfg.output_dir), "corrected")

        bh_cfg.__dict__["working_dir"] = bh_cfg.tmp_dir

        bh_cfg.__dict__["dataset"] = os.path.join(bh_cfg.output_dir,
            "dataset.info")

        if os.path.exists(bh_cfg.output_dir):
            shutil.rmtree(bh_cfg.output_dir)
        
        os.makedirs(bh_cfg.output_dir)
        if not os.path.exists(bh_cfg.working_dir):
            os.makedirs(bh_cfg.working_dir)

        log_filename = os.path.join(bh_cfg.output_dir, "correction.log")
        bh_cfg.__dict__["log_filename"] = log_filename

        print("\n===== Error correction started. Log can be found here: " +
              bh_cfg.log_filename + "\n")
        tee = support.Tee(log_filename, 'w', console=bh_cfg.output_to_console)

        if CONFIG_FILE:
            shutil.copy(CONFIG_FILE, bh_cfg.output_dir)

        # parsing dataset section
        bh_cfg.__dict__["single_cell"] = cfg["dataset"].single_cell
        bh_cfg.__dict__["paired_reads"] = []
        bh_cfg.__dict__["single_reads"] = []

        import bh_aux

        for k, v in cfg["dataset"].__dict__.items():
            if not isinstance(v, list):
                v = [v]

            # saving original reads to dataset
            if k.find("_reads") != -1:
                quoted_value = '"'
                for item in v:
                    quoted_value += os.path.abspath(os.path.expandvars(item)) + ' '
                quoted_value += '"'
                bh_cfg.__dict__["original_" + k] = quoted_value

            # saving reference to dataset in developer_mode
            if bh_cfg.developer_mode:
                if "reference" in cfg["dataset"].__dict__:
                    bh_cfg.__dict__["reference_genome"] = os.path.abspath(
                        os.path.expandvars(cfg["dataset"].reference))                       

            if k.startswith("single_reads"):
                for item in v:
                    item = os.path.abspath(os.path.expandvars(item))
                    item = bh_aux.ungzip_if_needed(item, bh_cfg.working_dir)
                    if not bh_cfg.single_reads:
                        bh_cfg.single_reads.append(item)
                    else:
                        bh_cfg.single_reads[0] = bh_aux.merge_single_files(item,
                            bh_cfg.single_reads[0], bh_cfg.working_dir)

            elif k.startswith("paired_reads"):
                cur_paired_reads = []
                if len(v) == 1:
                    item = os.path.abspath(os.path.expandvars(v[0]))
                    cur_paired_reads = bh_aux.split_paired_file(item, bh_cfg.working_dir)
                elif len(v) == 2:
                    for item in v:
                        item = os.path.abspath(os.path.expandvars(item))
                        item = bh_aux.ungzip_if_needed(item, bh_cfg.working_dir)
                        cur_paired_reads.append(item)

                if not bh_cfg.paired_reads:
                    bh_cfg.paired_reads = cur_paired_reads
                else:
                    bh_cfg.paired_reads = bh_aux.merge_paired_files(cur_paired_reads,
                        bh_cfg.paired_reads, bh_cfg.working_dir)

        bh_dataset_filename = bh_logic.run_bh(spades_home, execution_home, bh_cfg)

        tee.free()
        print("\n===== Error correction finished. Log can be found here: " +
              bh_cfg.log_filename + "\n")

    result_contigs_filename = ""
    result_scaffolds_filename = ""
    if "assembly" in cfg:
        spades_cfg = merge_configs(cfg["assembly"], cfg["common"])

        if "HEAPCHECK" in os.environ:
            del os.environ["HEAPCHECK"]
        if "heap_check" in spades_cfg.__dict__:
            os.environ["HEAPCHECK"] = spades_cfg.heap_check

        if "error_correction" in cfg:
            spades_cfg.__dict__["align_original_reads"] = True
        else:
            spades_cfg.__dict__["align_original_reads"] = False

        has_paired = False
        for k in cfg["dataset"].__dict__.keys():
            if k.startswith("paired_reads"):
                has_paired = True
                break
        if has_paired:
            spades_cfg.__dict__["paired_mode"] = True
        else:
            spades_cfg.__dict__["paired_mode"] = False

        spades_cfg.__dict__["log_filename"] = os.path.join(spades_cfg.output_dir,
            "assembly.log")
        spades_cfg.__dict__["result_contigs"] = os.path.join(spades_cfg.output_dir,
            "contigs.fasta")
        spades_cfg.__dict__["result_scaffolds"] = os.path.join(spades_cfg.output_dir,
            "scaffolds.fasta")

        spades_cfg.__dict__["additional_contigs"] = os.path.join(spades_cfg.output_dir,
            "simplified_contigs.fasta")

        print("\n===== Assembling started. Log can be found here: " + spades_cfg.log_filename +
              "\n")
        tee = support.Tee(spades_cfg.log_filename, 'w', console=spades_cfg.output_to_console)

        if CONFIG_FILE:
            shutil.copy(CONFIG_FILE, spades_cfg.output_dir)

        # dataset created during error correction
        if bh_dataset_filename:
            spades_cfg.__dict__["dataset"] = bh_dataset_filename

        if not "dataset" in spades_cfg.__dict__:
            # creating dataset
            dataset_filename = os.path.join(spades_cfg.output_dir, "dataset.info")
            dataset_file = open(dataset_filename, 'w')
            for k, v in cfg["dataset"].__dict__.iteritems():
                dataset_file.write(k + '\t')

                if isinstance(v, bool):
                    dataset_file.write(bool_to_str(v))
                else:
                    dataset_file.write('"')
                    if not isinstance(v, list):
                        v = [v]
                    for item in v:
                        item = os.path.abspath(os.path.expandvars(item))
                        dataset_file.write(str(item) + ' ')
                    dataset_file.write('"')
                dataset_file.write('\n')

            # saving reference to dataset in developer_mode
            if spades_cfg.developer_mode:
                if "reference" in cfg["dataset"].__dict__:
                    dataset_file.write("reference_genome" + '\t')
                    dataset_file.write(os.path.abspath(
                        os.path.expandvars(cfg["dataset"].reference)) + '\n')                    

            dataset_file.close()
            spades_cfg.__dict__["dataset"] = dataset_filename

        result_contigs_filename, result_scaffolds_filename = spades_logic.run_spades(spades_home, execution_home, spades_cfg)

        tee.free()
        print("\n===== Assembling finished. Log can be found here: " + spades_cfg.log_filename +
              "\n")

    print ("")
    if bh_dataset_filename:
        print (" * Corrected reads are in " + os.path.dirname(bh_dataset_filename) + "/")
    if result_contigs_filename:
        print (" * Assembled contigs are " + result_contigs_filename)
    if result_scaffolds_filename:
        print (" * Assembled scaffolds are " + result_scaffolds_filename)
    print ("")
    print ("Thank you for using SPAdes!")

    print("\n======= SPAdes pipeline finished\n")


if __name__ == '__main__':
    main()