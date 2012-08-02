#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import glob

import spades_init

spades_init.init()
spades_home = spades_init.spades_home

import support
from process_cfg import *

def error(err_str, prefix="== Error == "):
    print >> sys.stderr, "\n\n" + prefix + " " + err_str + "\n\n"
    exit(1)


def warning(warn_str, prefix="== Warning == "):
    print("\n\n" + prefix + " " + warn_str + "\n\n")


def prepare_config_bh(filename, cfg):
    subst_dict = dict()
    cfg.working_dir = os.path.abspath(cfg.working_dir)

    if len(cfg.paired_reads) == 2:
        subst_dict["input_paired_1"] = cfg.paired_reads[0]
        subst_dict["input_paired_2"] = cfg.paired_reads[1]
    if len(cfg.single_reads) == 1:
        subst_dict["input_single"] = cfg.single_reads[0]

    subst_dict["input_working_dir"] = cfg.working_dir
    subst_dict["general_max_iterations"] = cfg.max_iterations
    subst_dict["general_max_nthreads"] = cfg.max_threads
    subst_dict["count_merge_nthreads"] = cfg.max_threads
    subst_dict["bayes_nthreads"] = cfg.max_threads
    subst_dict["expand_nthreads"] = cfg.max_threads
    subst_dict["correct_nthreads"] = cfg.max_threads
    subst_dict["general_hard_memory_limit"] = cfg.max_memory

    if "qvoffset" in cfg.__dict__:
        subst_dict["input_qvoffset"] = cfg.qvoffset

    substitute_params(filename, subst_dict)


def prepare_config_spades(filename, cfg, prev_K, K, last_one):
    subst_dict = dict()
    cfg.working_dir = os.path.abspath(cfg.working_dir)

    subst_dict["K"] = str(K)
    subst_dict["dataset"] = cfg.dataset
    subst_dict["output_base"] = cfg.working_dir
    subst_dict["additional_contigs"] = cfg.additional_contigs
    subst_dict["entry_point"] = "construction"
    subst_dict["developer_mode"] = bool_to_str(cfg.developer_mode)
    subst_dict["SAM_writer_enable"] = bool_to_str(cfg.generate_sam_files and last_one)
    subst_dict["align_original_reads"] = bool_to_str(cfg.align_original_reads)
    subst_dict["align_before_RR"] = bool_to_str(not cfg.paired_mode)
    subst_dict["align_after_RR"] = bool_to_str(cfg.paired_mode)
    subst_dict["project_name"] = ""
    subst_dict["gap_closer_enable"] = bool_to_str(last_one and cfg.gap_closer)
    subst_dict["paired_mode"] = bool_to_str(last_one and cfg.paired_mode)
    subst_dict["additional_ec_removing"] = bool_to_str(last_one)
    subst_dict["use_additional_contigs"] = bool_to_str(prev_K)
    subst_dict["max_threads"] = cfg.max_threads
    subst_dict["max_memory"] = cfg.max_memory

    substitute_params(filename, subst_dict)


def print_used_values(cfg):
    def print_value(cfg, section, param, pretty_param="", margin="  "):
        if not pretty_param:
            pretty_param = param.capitalize().replace('_', ' ')
        sys.stdout.write(margin + pretty_param)
        if param in cfg[section].__dict__:
            print ": " + str(cfg[section].__dict__[param])
        else:
            if param.find("offset") != -1:
                print " will be auto-detected"

    print ""

    # main
    print_value(cfg, "common", "project_name", "", "")
    print_value(cfg, "common", "output_dir", "", "")
    print "Mode:",
    if ("error_correction" in cfg) and (not "assembly" in cfg):
        print "ONLY error correction (without assembler)"
    elif (not "error_correction" in cfg) and ("assembly" in cfg):
        print "ONLY assembler (without error correction)"
    else:
        print "error correction and assembler"
    if ("common" in cfg) and ("developer_mode" in cfg["common"].__dict__):
        if cfg["common"].developer_mode:
            print "Debug mode turned ON"
        else:
            print "Debug mode turned OFF"
    print ""

    # dataset
    if "dataset" in cfg:
        print "Dataset parameters:"

        if cfg["dataset"].single_cell:
            print "  Single-cell mode"
        else:
            print "  Multi-cell mode (you should set '--sc' flag if input data"\
                  " was obtained with MDA (single-cell) technology"

        no_single = True
        no_paired = True
        for k, v in cfg["dataset"].__dict__.iteritems():
            if k.startswith("single_reads") or k.startswith("paired_reads"):
                if k.startswith("paired_reads"):
                    no_paired = False
                    if not isinstance(v, list):
                        print "  Paired reads (file with interlaced reads):"
                    else:
                        print "  Paired reads (files with left and right reads):"
                else:
                    no_single = False
                    print "  Single reads:"
                if not isinstance(v, list):
                    v = [v]
                for reads_file in v:
                    print "    ", os.path.abspath(os.path.expandvars(reads_file))

        if no_paired:
            print "  Paired reads was not specified"
        if no_single:
            print "  Single reads was not specified"

    # error correction
    if "error_correction" in cfg:
        print "Error correction parameters:"
        print_value(cfg, "error_correction", "tmp_dir", "Dir for temp files")
        print_value(cfg, "error_correction", "max_iterations", "Iterations")
        print_value(cfg, "error_correction", "qvoffset", "PHRED offset")
        print "  Corrected reads will",
        if not cfg["error_correction"].gzip_output:
            print "NOT",
        print "be compressed (with gzip)"

    # assembly
    if "assembly" in cfg:
        print "Assembly parameters:"
        print_value(cfg, "assembly", "iterative_K", "k")
        print "  SAM file will",
        if not cfg["assembly"].generate_sam_files:
            print "NOT be generated (WARNING: SAM file are required "\
                  "for some of postprocessing tools)"
        else:
            print "be generated"
        print "  The gap closer will",
        if not cfg["assembly"].gap_closer:
            print "NOT",
        print "be used"

    print "Other parameters:"
    print_value(cfg, "common", "max_threads", "Threads")
    print_value(cfg, "common", "max_memory", "Memory limit (in Gb)", "  ")
    print ""


def check_config(cfg, default_project_name=""):
    ## special case: 0 iterations for the error correction means "No error correction!"

    if ("error_correction" in cfg) and ("max_iterations" in cfg["error_correction"].__dict__):
        if cfg["error_correction"].max_iterations == 0:
            del cfg["error_correction"]

    ## checking mandatory sections

    if (not "dataset" in cfg) and not (("assembly" in cfg) and ("dataset" in
                                                                cfg["assembly"].__dict__)):
        error("wrong config! You should specify 'dataset' section!")
        return False

    if (not "error_correction" in cfg) and (not "assembly" in cfg):
        error("wrong config! You should specify either 'error_correction' section (for reads"
              " error correction) or 'assembly' one (for assembling) or both!")
        return False

    ## checking existence of all files in dataset section

    no_files_with_reads = True
    for k, v in cfg["dataset"].__dict__.iteritems():
        if k.startswith("single_reads") or k.startswith("paired_reads"):
            no_files_with_reads = False
            if not isinstance(v, list):
                v = [v]
            for reads_file in v:
                if not os.path.isfile(os.path.expandvars(reads_file)):
                    error("file with reads doesn't exist! " + os.path.expandvars(reads_file))
                    return False
                else:
                    ext = os.path.splitext(os.path.expandvars(reads_file))[1]
                    if ext not in ['.fa', '.fasta', '.fq', '.fastq', '.gz']:
                        error("file with reads has unsupported format (only .fa, .fasta, .fq,"
                              " .fastq, .gz are supported)! " + os.path.expandvars(reads_file))
                        return False

    if no_files_with_reads:
        error("wrong config! You should specify at least one file with reads!")
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

    if not "project_name" in cfg["common"].__dict__:
        cfg["common"].__dict__["project_name"] = default_project_name

    cfg["common"].output_dir = os.path.join(os.path.abspath(
        os.path.expandvars(cfg["common"].output_dir)), cfg["common"].project_name)

    # dataset
    if "dataset" in cfg:
        if not "single_cell" in cfg["dataset"].__dict__:
            cfg["dataset"].__dict__["single_cell"] = False
        if cfg["common"].developer_mode and ("quality_assessment" in cfg) and ("reference" in cfg["quality_assessment"]):
            cfg["dataset"].__dict__["reference"] = cfg["quality_assessment"].reference

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

long_options = "12= threads= memory= tmp-dir= iterations= phred-offset= sc "\
               "generate-sam-file only-error-correction only-assembler "\
               "disable-gap-closer disable-gzip-output help test debug reference= help-hidden".split()
short_options = "n:o:1:2:s:k:t:m:i:h"


def check_file(f, message=''):
    if not os.path.isfile(f):
        error("file not found (%s): %s" % (message, f))
    return f


def usage(show_hidden=False):
    print >> sys.stderr, "SPAdes genome assembler"
    print >> sys.stderr, "Usage:", sys.argv[0], "[options] -n <project name>"
    print >> sys.stderr, ""
    print >> sys.stderr, "Options:"
    print >> sys.stderr, "-n\t<project_name>\tname of the project"
    print >> sys.stderr, "-o\t<output_dir>\tdirectory to store all result files"\
                         " [default: spades_output]"
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
    if show_hidden:
        print >> sys.stderr, "--reference\t<filename>\tfile with reference for deep analysis"\
                         " (only in debug mode)"
    print >> sys.stderr, "-t/--threads\t<int>\t\tnumber of threads [default: 16]"
    print >> sys.stderr, "-m/--memory\t<int>\t\tRAM limit for SPAdes in Gb"\
                         " (terminates if exceeded) [default: 250]"
    print >> sys.stderr, "--tmp-dir\t<dirname>\tdirectory for error correction's"\
                         " temp files"
    print >> sys.stderr, "\t\t\t\t[default: <output_dir>/<project_name>/corrected/tmp]"
    print >> sys.stderr, "-k\t<int,int,...>\t\tcomma-separated list of k-mer sizes"\
                         " (must be odd and less than 100) [default: 21,33,55]"
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

    print >> sys.stderr, ""
    print >> sys.stderr, "--test\t\trun SPAdes on toy dataset"
    print >> sys.stderr, "--debug\t\trun SPAdes in debug mode"
    print >> sys.stderr, "-h/--help\tprint this usage message"
    if show_hidden:
        print >> sys.stderr, "--help-hidden\tprint this usage message with all hidden options"

    print >> sys.stderr, ""
    print >> sys.stderr, "or you can run SPAdes with config file:", sys.argv[0],\
    "<config file>"
    print >> sys.stderr, "sample config is spades_config.info"


def main():
    os.environ["LC_ALL"] = "C"
    
    CONFIG_FILE = ""
    options = None

    if len(sys.argv) == 1:
        usage()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if os.path.isfile(sys.argv[1]):
            CONFIG_FILE = sys.argv[1]
    if not CONFIG_FILE:
        try:
            options, not_options = getopt.gnu_getopt(sys.argv, short_options, long_options)
        except getopt.GetoptError, err:
            print >> sys.stderr, err
            print >> sys.stderr
            usage()
            sys.exit(1)

    if not CONFIG_FILE and not options:
        usage()
        sys.exit(1)

    # all parameters are stored here
    cfg = dict()

    if options:
        project_name = ""
        output_dir = ""
        tmp_dir = ""
        reference = ""

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

        developer_mode = False

        threads = None
        memory = None
        qvoffset = None
        iterations = None

        for opt, arg in options:
            if opt == '-n':
                project_name = arg
            elif opt == '-o':
                output_dir = arg
            elif opt == "--tmp-dir":
                tmp_dir = arg
            elif opt == "--reference":
                reference = check_file(arg, 'reference')

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
                    if k > 100:
                        error('wrong k value ' + str(k) + ': all k values should be less than 100')
                    if k % 2 == 0:
                        error('wrong k value ' + str(k) + ': all k values should be odd')

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

            elif opt == '-h' or opt == "--help":
                usage()
                sys.exit(0)
            elif opt == "--help-hidden":
                usage(True)
                sys.exit(0)

            elif opt == "--test": # running test
                if os.path.isfile("spades_config.info"):
                    CONFIG_FILE = "spades_config.info"
                elif os.path.isfile(os.path.join(spades_home, "spades_config.info")):
                    CONFIG_FILE = os.path.join(spades_home, "spades_config.info")
                break
            else:
                raise ValueError

        if not CONFIG_FILE:
            if not project_name:
                error("the project name is not set! It is a mandatory parameter (-n PROJECT_NAME).")

            if len(paired1) != len(paired2):
                error("the number of files with left paired reads is not equal to the"
                      " number of files with right paired reads!")

            if not paired and not paired1 and not single:
                error("you should specify either paired reads (-1, -2 or -12) or single reads (-s) or both!")

            # filling cfg
            cfg["common"] = load_config_from_vars(dict())
            cfg["dataset"] = load_config_from_vars(dict())
            if not only_assembler:
                cfg["error_correction"] = load_config_from_vars(dict())
            if not only_error_correction:
                cfg["assembly"] = load_config_from_vars(dict())

            # filling reads
            paired_counter = 0
            if paired:
                for read in paired:
                    paired_counter += 1
                    cfg["dataset"].__dict__["paired_reads#" + str(paired_counter)] = read
                    #cfg["dataset"].__dict__["paired_reads"] = read

            if paired1:
                for i in range(len(paired1)):
                    paired_counter += 1
                    cfg["dataset"].__dict__["paired_reads#" + str(paired_counter)] = [paired1[i],
                                                                                      paired2[i]]
                    #cfg["dataset"].__dict__["paired_reads"] = [paired1[i], paired2[i]]  

            if single:
                cfg["dataset"].__dict__["single_reads"] = single
            
            # filling other parameters

            # common
            if output_dir:
                cfg["common"].__dict__["output_dir"] = output_dir
            if threads:
                cfg["common"].__dict__["max_threads"] = threads
            if memory:
                cfg["common"].__dict__["max_memory"] = memory
            if developer_mode:
                cfg["common"].__dict__["developer_mode"] = developer_mode

            # dataset
            cfg["dataset"].__dict__["single_cell"] = single_cell
            if developer_mode and reference:
                cfg["dataset"].__dict__["reference"] = reference

            # error correction
            if not only_assembler:
                if tmp_dir:
                    cfg["error_correction"].__dict__["tmp_dir"] = tmp_dir
                if qvoffset:
                    cfg["error_correction"].__dict__["qvoffset"] = qvoffset
                if iterations != None:
                    cfg["error_correction"].__dict__["max_iterations"] = iterations
                cfg["error_correction"].__dict__["gzip_output"] = not disable_gzip_output

            # assembly
            if not only_error_correction:
                if k_mers:
                    cfg["assembly"].__dict__["iterative_K"] = k_mers
                cfg["assembly"].__dict__["generate_sam_files"] = generate_sam_files
                cfg["assembly"].__dict__["gap_closer"] = not disable_gap_closer

            if not check_config(cfg, project_name):
                return

    if CONFIG_FILE:
        cfg = load_config_from_info_file(CONFIG_FILE)

        os.environ["cfg"] = os.path.dirname(os.path.abspath(CONFIG_FILE))

        if not check_config(cfg, os.path.splitext(os.path.basename(CONFIG_FILE))[0]):
            return

    print("\n======= SPAdes pipeline started\n")

    if not os.path.isdir(cfg["common"].output_dir):
        os.makedirs(cfg["common"].output_dir)
    log_filename = os.path.join(cfg["common"].output_dir, "params.txt")
    tee = support.Tee(log_filename, 'w', console=cfg["common"].output_to_console)

    if CONFIG_FILE:
        print("Using config file: " + CONFIG_FILE)
    else:
        print "Command line: ",
        for v in sys.argv:
            print v,
        print ""

    print_used_values(cfg)
    tee.free()

    bh_dataset_filename = ""
    if "error_correction" in cfg:
        bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])

        bh_cfg.output_dir = os.path.join(os.path.expandvars(bh_cfg.output_dir), "corrected")

        bh_cfg.__dict__["working_dir"] = bh_cfg.tmp_dir

        bh_cfg.__dict__["dataset"] = os.path.join(bh_cfg.output_dir,
            cfg["common"].project_name + ".dataset")

        start_bh = True
        if os.path.exists(bh_cfg.output_dir):
            if os.path.exists(bh_cfg.dataset):
                question = ["WARNING! It looks like error correction was already done!",
                            "Folder with corrected dataset " + bh_cfg.output_dir +
                            " already exists!",
                            "Do you want to overwrite this folder and start error"
                            " correction again?"]
                answer = support.question_with_timer(question, 10, 'n')
                if answer == 'n':
                    start_bh = False
                    bh_dataset_filename = bh_cfg.dataset
                    print("\n===== Error correction skipped\n")
                else:
                    os.remove(bh_cfg.dataset)
                    shutil.rmtree(bh_cfg.output_dir)

        if start_bh:
            if not os.path.exists(bh_cfg.output_dir):
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

            for k, v in cfg["dataset"].__dict__.iteritems():
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

            bh_dataset_filename = run_bh(bh_cfg)

            tee.free()
            print("\n===== Error correction finished. Log can be found here: " +
                  bh_cfg.log_filename + "\n")

    result_contigs_filename = ""
    if "assembly" in cfg:
        spades_cfg = merge_configs(cfg["assembly"], cfg["common"])

        if "error_correction" in cfg:
            spades_cfg.__dict__["align_original_reads"] = True
        else:
            spades_cfg.__dict__["align_original_reads"] = False

        has_paired = False
        for k in cfg["dataset"].__dict__.iterkeys():
            if k.startswith("paired_reads"):
                has_paired = True
                break
        if has_paired:
            spades_cfg.__dict__["paired_mode"] = True
        else:
            spades_cfg.__dict__["paired_mode"] = False

        def make_link(where, link):
            if os.path.islink(link):
                os.remove(link)
            if not os.path.exists(link):
                os.symlink(where, link)

        def make_working_dir(output_dir):
            import datetime

            name = "spades_" + datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
            working_dir = os.path.join(output_dir, name)
            os.makedirs(working_dir)
            return working_dir

        spades_cfg.__dict__["working_dir"] = make_working_dir(spades_cfg.output_dir)
        spades_cfg.__dict__["log_filename"] = os.path.join(spades_cfg.working_dir,
            "assembly.log")
        spades_cfg.__dict__["result_contigs"] = os.path.join(spades_cfg.working_dir,
            spades_cfg.project_name + ".fasta")
        spades_cfg.__dict__["additional_contigs"] = os.path.join(spades_cfg.working_dir,
            "simplified_contigs.fasta")
        final_contigs_folder = os.path.join(spades_cfg.output_dir, "contigs")

        make_link(os.path.basename(spades_cfg.working_dir), os.path.join(spades_cfg.output_dir,
            "latest"))
        if os.path.exists(final_contigs_folder):
            shutil.rmtree(final_contigs_folder)

        print("\n===== Assembling started. Log can be found here: " + spades_cfg.log_filename +
              "\n")
        tee = support.Tee(spades_cfg.log_filename, 'w', console=spades_cfg.output_to_console)

        if CONFIG_FILE:
            shutil.copy(CONFIG_FILE, spades_cfg.working_dir)

            # dataset created during error correction
        if bh_dataset_filename:
            spades_cfg.__dict__["dataset"] = bh_dataset_filename

        if not "dataset" in spades_cfg.__dict__:
            # creating dataset
            dataset_filename = os.path.join(spades_cfg.working_dir, cfg["common"].project_name +
                                                                    ".dataset")
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
        else:
            spades_cfg.dataset = os.path.abspath(os.path.expandvars(spades_cfg.dataset))
            shutil.copy(spades_cfg.dataset, spades_cfg.working_dir)
            # for developers: if dataset was set in 'assembly' section then use reference from it
            if ("quality_assessment" in cfg) and not bh_dataset_filename:
                for key_to_del in ["reference", "genes", "operons"]:
                    if key_to_del in cfg["quality_assessment"].__dict__:
                        del cfg["quality_assessment"].__dict__[key_to_del]

                dataset_cfg = load_config_from_file(spades_cfg.dataset)
                for k, v in dataset_cfg.__dict__.iteritems():
                    if k == "reference_genome":
                        cfg["quality_assessment"].__dict__["reference"] = os.path.join(
                            os.path.dirname(spades_cfg.dataset), v)

        result_contigs_filename = run_spades(spades_cfg)

        tee.free()
        print("\n===== Assembling finished. Log can be found here: " + spades_cfg.log_filename +
              "\n")

        make_link(os.path.basename(spades_cfg.working_dir), os.path.join(spades_cfg.output_dir,
            "latest_success"))
        if not os.path.exists(final_contigs_folder):
            os.makedirs(final_contigs_folder)
        shutil.copy(result_contigs_filename, final_contigs_folder)
        shutil.copy(spades_cfg.log_filename, final_contigs_folder)
        sam_file_linkname = os.path.splitext(result_contigs_filename)[0] + ".sam"
        if os.path.exists(sam_file_linkname):
            os.symlink(sam_file_linkname, os.path.join(final_contigs_folder,
                os.path.basename(sam_file_linkname)))

    quality_final_report = ""
    if ("quality_assessment" in cfg) and result_contigs_filename:
        quality_cfg = merge_configs(cfg["quality_assessment"], cfg["common"])

        quality_cfg.__dict__["working_dir"] = os.path.dirname(result_contigs_filename)
        quality_cfg.__dict__["log_filename"] = os.path.join(quality_cfg.working_dir,
            "quality.log")
        quality_cfg.__dict__["result_contigs"] = result_contigs_filename

        print("\n===== Quality assessment started. Log can be found here: " +
              quality_cfg.log_filename + "\n")
        tee = support.Tee(quality_cfg.log_filename, 'w', console=quality_cfg.output_to_console)

        quality_final_report = run_quality(quality_cfg)

        tee.free()
        print("\n===== Quality assessment finished. Log can be found here: " +
              quality_cfg.log_filename + "\n")

    print ""
    if bh_dataset_filename:
        print " * Corrected reads are in " + os.path.dirname(bh_dataset_filename) + "/"
    if result_contigs_filename:
        print " * Assembled contigs are " + result_contigs_filename
    if quality_final_report:
        print " * Assessment of their quality is in " + quality_final_report
    print ""
    print "Thank you for using SPAdes!"

    print("\n======= SPAdes pipeline finished\n")


def run_bh(cfg):
    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    shutil.copytree(os.path.join(spades_home, "configs"), dst_configs)
    cfg_file_name = os.path.join(dst_configs, "hammer", "config.info")

    prepare_config_bh(cfg_file_name, cfg)

    execution_home = os.path.join(spades_home, 'bin')
    command = os.path.join(execution_home,"hammer") + " " +\
              os.path.abspath(cfg_file_name)

    print("\n== Running error correction tool: " + command + "\n")
    support.sys_call(command)

    import bh_aux

    dataset_str = bh_aux.generate_dataset(cfg)
    dataset_filename = cfg.dataset
    dataset_file = open(dataset_filename, "w")
    dataset_file.write(dataset_str)
    dataset_file.close()
    print("\n== Dataset description file created: " + dataset_filename + "\n")

    shutil.rmtree(cfg.working_dir)

    return dataset_filename


def run_spades(cfg):
    if not isinstance(cfg.iterative_K, list):
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    execution_home = os.path.join(spades_home, 'bin')

    count = 0
    prev_K = None

    bin_reads_dir = os.path.join(cfg.working_dir, ".bin_reads")

    for K in cfg.iterative_K:
        count += 1

        dst_configs = os.path.join(cfg.working_dir, "config_K" + str(K))
        os.mkdir(dst_configs)
        dst_configs = os.path.join(dst_configs, "configs")
        shutil.copytree(os.path.join(spades_home, "configs"), dst_configs)
        cfg_file_name = os.path.join(dst_configs, "debruijn", "config.info")

        prepare_config_spades(cfg_file_name, cfg, prev_K, K,
            count == len(cfg.iterative_K))
        prev_K = K

        command = os.path.join(execution_home, "spades") + " " +\
                  os.path.abspath(cfg_file_name)

        if os.path.isdir(bin_reads_dir):
            if glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
                for cor_filename in glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
                    cor_index = cor_filename.rfind("_cor")
                    new_bin_filename = cor_filename[:cor_index] + cor_filename[cor_index + 4:]
                    shutil.move(cor_filename, new_bin_filename)

        print("\n== Running assembler: " + command + "\n")

        support.sys_call(command, execution_home)

        latest = os.path.join(cfg.working_dir, "K%d" % (K), "latest")
        latest = os.readlink(latest)
        latest = os.path.join("K%d" % (K), latest)
        os.symlink(latest, os.path.join(cfg.working_dir, "link_K%d" % (K)))
        latest = os.path.join(cfg.working_dir, latest)
        # python2.4 doesn't support os.path.relpath
        # os.symlink(os.path.relpath(latest, cfg.working_dir),
        #   os.path.join(cfg.working_dir, "link_K%d" % (K)))

    shutil.copyfile(os.path.join(latest, "final_contigs.fasta"), cfg.result_contigs)
    if cfg.developer_mode:
        # before repeat resolver contigs
        before_RR_contigs = os.path.join(os.path.dirname(cfg.result_contigs),
            "contigs_before_RR.fasta")
        shutil.copyfile(os.path.join(latest, "contigs_before_RR.fasta"), before_RR_contigs)
        # saves
        os.symlink(os.path.join(latest, "saves"), os.path.join(
            os.path.dirname(cfg.result_contigs), "saves"))

    os.remove(cfg.additional_contigs)

    if glob.glob(os.path.join(latest, "*.sam")):
        sam_file_linkname = os.path.join(os.path.dirname(cfg.result_contigs),
            cfg.project_name + ".sam")
        os.symlink(glob.glob(os.path.join(latest, "*.sam"))[0], sam_file_linkname)

    if os.path.isdir(bin_reads_dir):
        shutil.rmtree(bin_reads_dir)

    return cfg.result_contigs


def run_quality(cfg):        

    args = [cfg.result_contigs]

    if cfg.developer_mode:
        before_RR_contigs = os.path.join(os.path.dirname(cfg.result_contigs),
            "contigs_before_RR.fasta")
        args.append(before_RR_contigs)

    if "reference" in cfg.__dict__:
        args.append("-R")
        args.append(os.path.abspath(os.path.expandvars(cfg.reference)))
    if "genes" in cfg.__dict__:
        args.append("-G")
        args.append(os.path.abspath(os.path.expandvars(cfg.genes)))
    if "operons" in cfg.__dict__:
        args.append("-O")
        args.append(os.path.abspath(os.path.expandvars(cfg.operons)))
    quality_output_dir = os.path.join(cfg.working_dir, "quality_results")
    args.append("-o")
    args.append(quality_output_dir)
    import quality

    quality.main(args, lib_dir=os.path.join(spades_home, "src/tools/quality/libs"))

    return os.path.join(quality_output_dir, "report.txt")


if __name__ == '__main__':
    try:
        main()
    except support.spades_error, e:
        print(e.what())
