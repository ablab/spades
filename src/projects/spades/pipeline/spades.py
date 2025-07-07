#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
import shutil
import platform
import sys
from site import addsitedir

import spades_init

spades_init.init()
spades_home = spades_init.spades_home
bin_home = spades_init.bin_home
python_modules_home = spades_init.python_modules_home
ext_python_modules_home = spades_init.ext_python_modules_home
spades_version = spades_init.spades_version

import support

support.check_python_version()

addsitedir(ext_python_modules_home)
import pyyaml3 as pyyaml

import options_storage
options_storage.spades_version = spades_version

import options_parser
from stages.pipeline import Pipeline
import executor_save_yaml


log = logging.getLogger("spades")


def print_used_values(cfg):
    def print_value(config, section, param, pretty_param="", margin="  "):
        if not pretty_param:
            pretty_param = param.capitalize().replace('_', ' ')
        line = margin + pretty_param
        if param in config[section].__dict__:
            line += ": " + str(config[section].__dict__[param])
        else:
            if "offset" in param:
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
    if not options_storage.args.rna:
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

        if options_storage.args.iontorrent:
            log.info("  IonTorrent data")
        if options_storage.args.sewage:
            log.info("  Sewage mode")
        if options_storage.args.bio:
            log.info("  BiosyntheticSPAdes mode")
        if options_storage.args.rnaviral:
            log.info("  RNA virus assembly mode")
        if options_storage.args.meta:
            log.info("  Metagenomic mode")
        elif options_storage.args.large_genome:
            log.info("  Large genome mode")
        elif options_storage.args.isolate:
            log.info("  Isolate mode")
        elif options_storage.args.rna:
            log.info("  RNA-seq mode")
        elif options_storage.args.single_cell:
            log.info("  Single-cell mode")
        else:
            log.info("  Standard mode")
            log.info("  For multi-cell/isolate data we recommend to use '--isolate' option;"
                     " for single-cell MDA data use '--sc';"
                     " for metagenomic data use '--meta';"
                     " for RNA-Seq use '--rna'.")

        log.info("  Reads:")
        dataset_data = pyyaml.load(open(cfg["dataset"].yaml_filename))
        dataset_data = support.relative2abs_paths(dataset_data, os.path.dirname(cfg["dataset"].yaml_filename))
        support.pretty_print_reads(dataset_data)

    # error correction
    if "error_correction" in cfg:
        log.info("Read error correction parameters:")
        print_value(cfg, "error_correction", "max_iterations", "Iterations")
        print_value(cfg, "error_correction", "qvoffset", "PHRED offset")

        if cfg["error_correction"].gzip_output:
            log.info("  Corrected reads will be compressed")
        else:
            log.info("  Corrected reads will NOT be compressed")

    # assembly
    if "assembly" in cfg:
        log.info("Assembly parameters:")
        if options_storage.auto_K_allowed():
            log.info("  k: automatic selection based on read length")
        else:
            print_value(cfg, "assembly", "iterative_K", "k")
        if options_storage.args.plasmid:
            log.info("  Extrachromosomal mode is turned ON")
        if cfg["assembly"].disable_rr:
            log.info("  Repeat resolution is DISABLED")
        else:
            log.info("  Repeat resolution is enabled")
        if options_storage.args.careful:
            log.info("  Mismatch careful mode is turned ON")
        else:
            log.info("  Mismatch careful mode is turned OFF")
        if "mismatch_corrector" in cfg:
            log.info("  MismatchCorrector will be used")
        else:
            log.info("  MismatchCorrector will be SKIPPED")
        if cfg["assembly"].cov_cutoff == "off":
            log.info("  Coverage cutoff is turned OFF")
        elif cfg["assembly"].cov_cutoff == "auto":
            log.info("  Coverage cutoff is turned ON and threshold will be auto-detected")
        else:
            log.info("  Coverage cutoff is turned ON and threshold is %f" % cfg["assembly"].cov_cutoff)
        if cfg["assembly"].gfa11:
            log.info("  Assembly graph output will use GFA v1.1 format")
        else:
            log.info("  Assembly graph output will use GFA v1.2 format")

    log.info("Other parameters:")
    print_value(cfg, "common", "tmp_dir", "Dir for temp files")
    print_value(cfg, "common", "max_threads", "Threads")
    print_value(cfg, "common", "max_memory", "Memory limit (in Gb)", "  ")
    log.info("")


def create_logger():
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)


def check_cfg_for_partial_run(cfg, partial_run_type="restart-from"):  # restart-from ot stop-after
    if partial_run_type == "restart-from":
        check_point = options_storage.args.restart_from
        action = "restart from"
        verb = "was"
    elif partial_run_type == "stop-after":
        check_point = options_storage.args.stop_after
        action = "stop after"
        verb = "is"
    else:
        return

    if check_point == "ec" and ("error_correction" not in cfg):
        support.error(
            "failed to %s 'read error correction' ('%s') because this stage %s not specified!" % (action, check_point, verb))
    if check_point == "mc" and ("mismatch_corrector" not in cfg):
        support.error(
            "failed to %s 'mismatch correction' ('%s') because this stage %s not specified!" % (action, check_point, verb))
    if check_point == "as" or check_point.startswith('k'):
        if "assembly" not in cfg:
            support.error(
                "failed to %s 'assembling' ('%s') because this stage %s not specified!" % (action, check_point, verb))

def get_options_from_params(params_filename, running_script):
    command_line = None
    options = None
    prev_running_script = None
    if not os.path.isfile(params_filename):
        return command_line, options, prev_running_script, \
               "failed to parse command line of the previous run (%s not found)!" % params_filename

    with open(params_filename) as params:
        command_line = params.readline().strip()
        spades_prev_version = None
        for line in params:
            if "SPAdes version:" in line:
                spades_prev_version = line.split("SPAdes version:")[1]
                break

    if spades_prev_version is None:
        return command_line, options, prev_running_script, \
               "failed to parse SPAdes version of the previous run!"
    if spades_prev_version.strip() != spades_version.strip():
        return command_line, options, prev_running_script, \
               "SPAdes version of the previous run (%s) is not equal to the current version of SPAdes (%s)!" \
               % (spades_prev_version.strip(), spades_version.strip())
    if "Command line: " not in command_line or '\t' not in command_line:
        return command_line, options, prev_running_script, "failed to parse executable script of the previous run!"
    options = command_line.split('\t')[1:]
    prev_running_script = command_line.split('\t')[0][len("Command line: "):]
    prev_running_script = os.path.basename(prev_running_script)
    running_script = os.path.basename(running_script)
    # we cannot restart/continue spades.py run with metaspades.py/rnaspades.py/etc and vice versa
    if prev_running_script != running_script:
        message = "executable script of the previous run (%s) is not equal " \
                  "to the current executable script (%s)!" % (prev_running_script, running_script)
        return command_line, options, prev_running_script, message
    return command_line, options, prev_running_script, ""


# parse options and safe all parameters to cfg
def parse_args(args):
    options, cfg, dataset_data = options_parser.parse_args(bin_home, spades_home,
                                                           secondary_filling=False, restart_from=False)

    command_line = ""

    if options_storage.args.continue_mode:
        restart_from = options_storage.args.restart_from
        command_line, options, script_name, err_msg = get_options_from_params(
            os.path.join(options_storage.args.output_dir, "params.txt"),
            args[0])
        if err_msg:
            support.error(err_msg + " Please restart from the beginning or specify another output directory.")
        options, cfg, dataset_data = options_parser.parse_args(bin_home, spades_home, secondary_filling=True,
                                                               restart_from=(options_storage.args.restart_from is not None),
                                                               options=options)

        options_storage.args.continue_mode = True
        options_storage.args.restart_from = restart_from

        if options_storage.args.restart_from:
            check_cfg_for_partial_run(cfg, partial_run_type="restart-from")

    if options_storage.args.stop_after:
        check_cfg_for_partial_run(cfg, partial_run_type="stop-after")

    support.check_single_reads_in_options()
    return cfg, dataset_data, command_line


def add_file_to_log(cfg, log_intance):
    log_filename = os.path.join(cfg["common"].output_dir, "spades.log")
    if options_storage.args.continue_mode:
        log_handler = logging.FileHandler(log_filename, mode='a')
    else:
        log_handler = logging.FileHandler(log_filename, mode='w')
    log_intance.addHandler(log_handler)
    return log_filename, log_handler


def get_restart_from_command_line(args):
    updated_params = ""
    for i in range(1, len(args)):
        if not args[i].startswith("-o") and not args[i].startswith("--restart-from") and \
                        args[i - 1] != "-o" and args[i - 1] != "--restart-from":
            updated_params += "\t" + args[i]

    updated_params = updated_params.strip()
    restart_from_update_message = "Restart-from=" + options_storage.args.restart_from + "\n"
    restart_from_update_message += "with updated parameters: " + updated_params
    return updated_params, restart_from_update_message


def get_command_line(args):
    command = ""
    for v in args:
        # substituting relative paths with absolute ones (read paths, output dir path, etc)
        v, prefix = support.get_option_prefix(v)
        if v in options_storage.dict_of_rel2abs.keys():
            v = options_storage.dict_of_rel2abs[v]
        if prefix:
            command += prefix + ":"
        command += v + "\t"
    return command


def print_params(log_filename, command_line, args, cfg):
    if options_storage.args.continue_mode:
        log.info("\n======= SPAdes pipeline continued. Log can be found here: " + log_filename + "\n")
        log.info("Restored from " + command_line)
        log.info("")

    params_filename = os.path.join(cfg["common"].output_dir, "params.txt")
    params_handler = logging.FileHandler(params_filename, mode='w')
    log.addHandler(params_handler)

    if not options_storage.args.continue_mode:
        log.info("Command line: " + get_command_line(args))
    elif options_storage.args.restart_from:
        update_params, restart_from_update_message = get_restart_from_command_line(args)
        command_line += "\t" + update_params
        log.info(command_line)
        log.info(restart_from_update_message)
    else:
        log.info(command_line)

    print_used_values(cfg)
    log.removeHandler(params_handler)


def clear_configs(cfg, command_before_restart_from, stage_id_before_restart_from):
    def matches_with_restart_from_arg(stage, restart_from_arg):
        return stage["short_name"].startswith(restart_from_arg.split(":")[0])

    spades_commands_fpath = os.path.join(cfg["common"].output_dir, "run_spades.yaml")
    with open(spades_commands_fpath) as stream:
        old_pipeline = pyyaml.load(stream)

    restart_from_stage_id = None
    for num in range(len(old_pipeline)):
        stage = old_pipeline[num]
        if matches_with_restart_from_arg(stage, options_storage.args.restart_from):
            restart_from_stage_id = num
            break

    if (command_before_restart_from is not None and
            old_pipeline[stage_id_before_restart_from]["short_name"] != command_before_restart_from.short_name):
        support.error("new and old pipelines have difference before %s" % options_storage.args.restart_from, log)

    if command_before_restart_from is None:
        first_del = 0
    else:
        first_del = stage_id_before_restart_from + 1

    if restart_from_stage_id is not None:
        stage_filename = options_storage.get_stage_filename(restart_from_stage_id, old_pipeline[restart_from_stage_id]["short_name"])
        if os.path.isfile(stage_filename):
            os.remove(stage_filename)

    for delete_id in range(first_del, len(old_pipeline)):
        stage_filename = options_storage.get_stage_filename(delete_id, old_pipeline[delete_id]["short_name"])
        if os.path.isfile(stage_filename):
            os.remove(stage_filename)

        cfg_dir = old_pipeline[delete_id]["config_dir"]
        if cfg_dir != "" and os.path.isdir(os.path.join(cfg["common"].output_dir, cfg_dir)):
            shutil.rmtree(os.path.join(cfg["common"].output_dir, cfg_dir))


def get_first_incomplete_command(filename):
    with open(filename) as stream:
        old_pipeline = pyyaml.load(stream)

    first_incomplete_stage_id = 0
    while first_incomplete_stage_id < len(old_pipeline):
        stage_filename = options_storage.get_stage_filename(first_incomplete_stage_id, old_pipeline[first_incomplete_stage_id]["short_name"])
        if not os.path.isfile(stage_filename):
            return old_pipeline[first_incomplete_stage_id]
        first_incomplete_stage_id += 1

    return None


def get_command_and_stage_id_before_restart_from(draft_commands):
    restart_from_stage_name = options_storage.args.restart_from.split(":")[0]

    if options_storage.args.restart_from == options_storage.LAST_STAGE:
        last_command = get_first_incomplete_command(os.path.join(get_stage.cfg["common"].output_dir, "run_spades.yaml"))
        if last_command is None:
            restart_from_stage_name = draft_commands[-1].short_name
        else:
            restart_from_stage_name = last_command["short_name"]

    restart_from_stage_id = None
    for num in range(len(draft_commands)):
        stage = draft_commands[num]
        if stage.short_name.startswith(restart_from_stage_name):
            restart_from_stage_id = num
            break

    if restart_from_stage_id is None:
        support.error(
            "failed to restart from %s because this stage was not specified!" % options_storage.args.restart_from,
            log)

    if ":" in options_storage.args.restart_from or options_storage.args.restart_from == options_storage.LAST_STAGE:
        return draft_commands[restart_from_stage_id], restart_from_stage_id

    if restart_from_stage_id > 0:
        stage_filename = options_storage.get_stage_filename(restart_from_stage_id - 1, draft_commands[restart_from_stage_id - 1].short_name)
        if not os.path.isfile(stage_filename):
            support.error(
                "cannot restart from stage %s: previous stage was not complete." % options_storage.args.restart_from,
                log)
        return draft_commands[restart_from_stage_id - 1], restart_from_stage_id - 1
    return None, None


def print_info_about_output_files(cfg, output_files):
    def check_and_report_output_file(output_file_key, message_prefix_text, error_message = ""):
        if os.path.isfile(output_files[output_file_key]):
            message = message_prefix_text + support.process_spaces(output_files[output_file_key])
            log.info(message)
        else:
            if error_message != "":
                log.info(error_message)

    if "error_correction" in cfg and os.path.isdir(
            os.path.dirname(output_files["corrected_dataset_yaml_filename"])):
        log.info(" * Corrected reads are in " + support.process_spaces(
            os.path.dirname(output_files["corrected_dataset_yaml_filename"]) + "/"))

    if "assembly" in cfg:
        error_message = ""
        if options_storage.args.plasmid:
            error_message = "No plasmid contigs assembled!!"
            if options_storage.args.meta:
                error_message = "No complete extrachromosomal contigs assembled!!"
        check_and_report_output_file("result_contigs_filename", " * Assembled contigs are in ", error_message)

        if options_storage.args.bio or options_storage.args.custom_hmms or options_storage.args.corona:
            check_and_report_output_file("result_domain_graph_filename", " * Domain graph is in ")
            check_and_report_output_file("result_gene_clusters_filename", " * Gene cluster sequences are in ")
            check_and_report_output_file("result_bgc_stats_filename", " * BGC cluster statistics ")

        if options_storage.args.rna:
            check_and_report_output_file("result_transcripts_filename", " * Assembled transcripts are in ")
            check_and_report_output_file("result_transcripts_paths_filename",
                                         " * Paths in the assembly graph corresponding to the transcripts are in ")

            for filtering_type in options_storage.filtering_types:
                result_filtered_transcripts_filename = os.path.join(cfg["common"].output_dir,
                                                                    filtering_type + "_filtered_" +
                                                                    options_storage.transcripts_name)
                if os.path.isfile(result_filtered_transcripts_filename):
                    message = " * " + filtering_type.capitalize() + " filtered transcripts are in " + \
                              support.process_spaces(result_filtered_transcripts_filename)
                    log.info(message)
        else:
            check_and_report_output_file("result_scaffolds_filename", " * Assembled scaffolds are in ")
            check_and_report_output_file("result_contigs_paths_filename",
                                         " * Paths in the assembly graph corresponding to the contigs are in ")
            check_and_report_output_file("result_scaffolds_paths_filename",
                                         " * Paths in the assembly graph corresponding to the scaffolds are in ")

        check_and_report_output_file("result_assembly_graph_filename", " * Assembly graph is in ")
        check_and_report_output_file("result_assembly_graph_filename_gfa", " * Assembly graph in GFA format is in ")


def get_output_files(cfg):
    output_files = dict()
    output_files["corrected_dataset_yaml_filename"] = ""
    output_files["result_contigs_filename"] = os.path.join(cfg["common"].output_dir, options_storage.contigs_name)
    output_files["result_scaffolds_filename"] = os.path.join(cfg["common"].output_dir, options_storage.scaffolds_name)

    output_files["result_assembly_graph_filename"] = os.path.join(cfg["common"].output_dir,
                                                                  options_storage.assembly_graph_name)
    output_files["result_assembly_graph_filename_gfa"] = os.path.join(cfg["common"].output_dir,
                                                                      options_storage.assembly_graph_name_gfa)
    output_files["result_contigs_paths_filename"] = os.path.join(cfg["common"].output_dir,
                                                                 options_storage.contigs_paths)
    output_files["result_scaffolds_paths_filename"] = os.path.join(cfg["common"].output_dir,
                                                                   options_storage.scaffolds_paths)
    output_files["result_transcripts_filename"] = os.path.join(cfg["common"].output_dir,
                                                               options_storage.transcripts_name)
    output_files["result_transcripts_paths_filename"] = os.path.join(cfg["common"].output_dir,
                                                                     options_storage.transcripts_paths)
    output_files["result_bgc_stats_filename"] = os.path.join(cfg["common"].output_dir, options_storage.bgc_stats_name)
    output_files["result_domain_graph_filename"] = os.path.join(cfg["common"].output_dir, options_storage.domain_graph_name)
    output_files["result_sewage_lineages_filename"] = os.path.join(cfg["common"].output_dir, options_storage.sewage_lineages)
    output_files["result_gene_clusters_filename"] = os.path.join(cfg["common"].output_dir, options_storage.scaffolds_name)
    output_files["misc_dir"] = os.path.join(cfg["common"].output_dir, "misc")
    # if mismatch correction is enabled then result contigs are copied to misc directory
    output_files["assembled_contigs_filename"] = os.path.join(output_files["misc_dir"], "assembled_contigs.fasta")
    output_files["assembled_scaffolds_filename"] = os.path.join(output_files["misc_dir"], "assembled_scaffolds.fasta")
    if options_storage.hmm_mode():
        output_files["result_scaffolds_filename"] = os.path.join(cfg["common"].output_dir, options_storage.secondary_scaffolds_name)
        output_files["result_scaffolds_paths_filename"] = os.path.join(cfg["common"].output_dir,
                                                                       options_storage.secondary_scaffolds_paths)
        output_files["result_contigs_filename"] = os.path.join(cfg["common"].output_dir,
                                                               options_storage.secondary_contigs_name)

    return output_files


def get_stage(iteration_name):
    if not options_storage.args.continue_mode:
        return options_storage.BASE_STAGE

    if (options_storage.args.restart_from is not None and
            options_storage.args.restart_from != options_storage.LAST_STAGE):
        if (":" in options_storage.args.restart_from and
                iteration_name == options_storage.args.restart_from.split(":")[0]):
            return options_storage.args.restart_from.split(":")[-1]
        else:
            return options_storage.BASE_STAGE

    if get_stage.restart_stage is None:
        last_command = get_first_incomplete_command(os.path.join(get_stage.cfg["common"].output_dir, "run_spades.yaml"))

        if last_command is not None:
            get_stage.restart_stage = last_command["short_name"]
        else:
            get_stage.restart_stage = "finish"

    if iteration_name == get_stage.restart_stage:
        return options_storage.LAST_STAGE
    else:
        return options_storage.BASE_STAGE


def build_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                   ext_python_modules_home, python_modules_home):
    from stages import before_start_stage
    from stages import error_correction_stage
    from stages import spades_stage
    from stages import correction_stage
    from stages import check_test_stage
    from stages import breaking_scaffolds_stage
    from stages import preprocess_reads_stage
    from stages import terminating_stage

    before_start_stage.add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                       ext_python_modules_home, python_modules_home)
    preprocess_reads_stage.add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                           ext_python_modules_home, python_modules_home)
    error_correction_stage.add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                           ext_python_modules_home, python_modules_home)

    get_stage.cfg, get_stage.restart_stage = cfg, None
    spades_stage.add_to_pipeline(pipeline, get_stage, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                 ext_python_modules_home, python_modules_home)
    correction_stage.add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                     ext_python_modules_home, python_modules_home)
    check_test_stage.add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                     ext_python_modules_home, python_modules_home)
    breaking_scaffolds_stage.add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                             ext_python_modules_home, python_modules_home)
    terminating_stage.add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, bin_home,
                                         ext_python_modules_home, python_modules_home)

def get_executor():
    import importlib
    module_name = "executor_" + options_storage.args.grid_engine
    executor_module = importlib.import_module(module_name)
    return executor_module.Executor(log)


def get_sh_dump_executor():
    if options_storage.args.grid_engine == "local":
        return executor_save_yaml.Executor(log)
    else:
        import executor_save_mpi_sh
        return executor_save_mpi_sh.Executor(log)


def check_dir_is_empty(dir_name):
    if dir_name is not None and \
            os.path.exists(dir_name) and \
            os.listdir(dir_name):
        support.warning("output dir is not empty! Please, clean output directory before run.")


def init_parser(args):
    if options_parser.is_first_run():
        options_storage.first_command_line = args
        check_dir_is_empty(options_parser.get_output_dir_from_args())
    else:
        output_dir = options_parser.get_output_dir_from_args()
        if output_dir is None:
            support.error("the output_dir is not set! It is a mandatory parameter (-o output_dir).")

        command_line, options, script, err_msg = get_options_from_params(
            os.path.join(output_dir, "params.txt"),
            args[0])

        if err_msg != "":
            support.error(err_msg)

        options_storage.first_command_line = [script] + options


def main(args):
    os.environ["LC_ALL"] = "C"

    init_parser(args)

    if len(args) == 1:
        options_parser.usage(spades_version)
        sys.exit(0)

    jobs = []
    executor = None
    pipeline = Pipeline()

    create_logger()
    cfg, dataset_data, command_line = parse_args(args)
    log_filename, log_handler = add_file_to_log(cfg, log)
    print_params(log_filename, command_line, args, cfg)

    if not options_storage.args.continue_mode:
        log.info("\n======= SPAdes pipeline started. Log can be found here: " + log_filename + "\n")

    support.check_binaries(bin_home)
    try:
        output_files = get_output_files(cfg)
        tmp_configs_dir = os.path.join(cfg["common"].output_dir, "configs")

        build_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data,
                       bin_home, ext_python_modules_home, python_modules_home)

        if options_storage.args.restart_from:
            draft_commands = pipeline.get_commands(cfg)
            command_before_restart_from, stage_id_before_restart_from = \
                get_command_and_stage_id_before_restart_from(draft_commands)
            clear_configs(cfg, command_before_restart_from, stage_id_before_restart_from)

        pipeline.generate_configs(cfg, spades_home, tmp_configs_dir)
        commands = pipeline.get_commands(cfg)

        executor = get_sh_dump_executor()
        executor.execute(commands)

        executor = get_executor()
        if options_storage.args.grid_engine != "local":
            executor.dump_commands(commands, os.path.join(options_storage.args.output_dir, "run_spades_on_cluster.sh"))

        if options_storage.args.only_generate_config:
            jobs = None
        else:
            jobs = executor.execute(commands)

        if jobs is not None and len(jobs):
            last_job = jobs[-1]
            if options_storage.args.grid_wait:
                log.info("Waiting for the last job: " + last_job)
                executor.join(last_job)
            else:
                log.info("Last job name: " + last_job)

        is_result = (options_storage.args.grid_engine != "save_yaml" and (
            not options_storage.args.only_generate_config)) and \
                    (options_storage.args.grid_wait or options_storage.args.grid_engine == "local")
        if is_result:
            # TODO make it executor method executor.is_fake()
            print_info_about_output_files(cfg, output_files)

        if not support.log_warnings(log):
            if is_result:
                log.info("\n======= SPAdes pipeline finished.")  # otherwise it finished WITH WARNINGS
            else:
                log.info("\n======= SPAdes pipeline submitted.")

    except KeyboardInterrupt:
        log.info("Ctrl + C pressed, killing jobs...")
        if executor is not None:
            for job in jobs:
                log.info("Killing " + job)
                executor.kill(job)
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            import errno
            if exc_type == OSError and exc_value.errno == errno.ENOEXEC:  # Exec format error
                support.error("it looks like you are using SPAdes binaries for another platform.\n" +
                              support.get_spades_binaries_info_message())
            else:
                log.exception(exc_value)
                support.error("exception caught: %s" % exc_type, log)
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            support.error("exception caught: %s" % exc_type, log)
    finally:
        mode = options_parser.get_mode()
        log.info("\nSPAdes log can be found here: %s" % log_filename)
        log.info("")

        if mode == "plasmid":
            spades_name = "plasmidSPAdes"
            cite_message = ("Antipov, D., Hartwick, N., Shen, M., Raiko, M., Lapidus, A. and Pevzner, P.A., 2016. "
                            "plasmidSPAdes: assembling plasmids from whole genome sequencing data. "
                            "Bioinformatics, 32(22), pp.3380-3387.")
            doi = "doi.org/10.1093/bioinformatics/btw493"
        elif mode == "bgc":
            spades_name = "BiosyntheticSPAdes"
            cite_message = ("Meleshko, D., Mohimani, H., Tracanna, V., Hajirasouliha, I., Medema, M.H., "
                            "Korobeynikov, A.  and Pevzner, P.A., 2019. BiosyntheticSPAdes: reconstructing "
                            "biosynthetic gene clusters from assembly graphs. Genome research, 29(8), pp.1352-1362.")
            doi = "doi.org/10.1101/gr.243477.118"
        elif mode == "meta":
            spades_name = "metaSPAdes"
            cite_message = ("Nurk, S., Meleshko, D., Korobeynikov, A. and Pevzner, P.A., 2017. "
                            "metaSPAdes: a new versatile metagenomic assembler. "
                            "Genome research, 27(5), pp.824-834.")
            doi = "doi.org/10.1101/gr.213959.116"
        elif mode == "metaplasmid":
            spades_name = "metaplasmidSPAdes"
            cite_message = ("Antipov, D., Raiko, M., Lapidus, A. and Pevzner, P.A., 2019. "
                            "Plasmid detection and assembly in genomic and metagenomic data sets. "
                            "Genome research, 29(6), pp.961-968.")
            doi = "doi.org/10.1101/gr.241299.118"
        elif mode == "metaviral":
            spades_name = "metaviralSPAdes"
            cite_message = ("Antipov, D., Raiko, M., Lapidus, A. and Pevzner, P.A., 2020. "
                            "Metaviral SPAdes: assembly of viruses from metagenomic data. "
                            "Bioinformatics, 36(14), pp.4126-4129.")
            doi = "doi.org/10.1093/bioinformatics/btaa490"
        elif mode == "corona":
            spades_name = "coronaSPAdes"
            cite_message = ("Meleshko, D., Hajirasouliha, I. and Korobeynikov, A., 2022. "
                            "coronaSPAdes: from biosynthetic gene clusters to RNA viral assemblies. "
                            "Bioinformatics, 38(1), pp.1-8.")
            doi = "doi.org/10.1093/bioinformatics/btab597"
        elif mode == "rnaviral":
            spades_name = "rnaviralSPAdes"
            cite_message = ("Korobeynikov, A. and Meleshko, D., 2022. Benchmarking state-of-the-art approaches "
                            "for norovirus genome assembly in metagenome sample.")
            doi = "doi.org/10.21203/rs.3.rs-1827448/v1"
        elif mode == "rna":
            spades_name = "rnaSPAdes"
            cite_message = ("Bushmanova, E., Antipov, D., Lapidus, A. and Prjibelski, A.D., 2019. "
                            "rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. "
                            "GigaScience, 8(9), p.giz100.")
            doi = "doi.org/10.1093/gigascience/giz100"
        else:
            spades_name = "SPAdes"
            cite_message = ("Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A. and Korobeynikov, A., 2020. "
                            "Using SPAdes de novo assembler. Current protocols in bioinformatics, 70(1), p.e102.")
            doi = "doi.org/10.1002/cpbi.102"

        log.info("Thank you for using %s! If you use it in your research, please cite:" % spades_name)
        log.info("")
        log.info("  %s" % cite_message)
        log.info("  %s" % doi)
        log.info("")

        log.removeHandler(log_handler)


if __name__ == "__main__":
    main(sys.argv)
