#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#script for testing SPAdes
#provide a path to .info file


import sys
import os
import shutil
import getopt
import re
import datetime
import argparse
import subprocess
import glob
import math
from traceback import print_exc

sys.path.append('./src/spades_pipeline/')
import process_cfg

sys.path.append('./src/test/teamcity/')
import teamcity_support
from teamcity_support import *

def parse_args_ut():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", "-c", help="teamcity.py info config OR", type=str)
    parser.set_defaults(config="")
    parser.add_argument("--config_dir", "-d", help="folder with configs, will be scanned recursively; not compatible with etalon contigs options", type=str)
    parser.set_defaults(config_dir="")
    parser.add_argument("--tmp_dir", "-t", help="folder for temporary QUAST results, will be cleared before and removed after the run (default ./.tc_tmp)", type=str)
    parser.set_defaults(tmp_dir="./.tc_tmp")
    parser.add_argument("--preserve_tmp_dir", "-p", help="do not clean or remove temporary dir", action='store_true')
    parser.set_defaults(preserve_tmp_dir=False)
    #all - update all; add all missibg metrics
    #basic - update only basic metrics; add only basic missing metrics
    #none - update only metrics that are already on config, add none
    parser.add_argument("--add_metrics", "-a", help="add new metrics to the info file [all | basic | none] (default: basic). none: do not add new thresholds to config; basic: add new thresholds only if they are from list of basic ones; all: add all thesholds", type=str)
    parser.set_defaults(add_metrics="basic")
    parser.add_argument("--etalon_contigs_prefix", help="prefix to saved contigs that will be used to compute thresholds; latest contigs from the config storage by default", type=str)
    parser.add_argument("--etalon_contigs_suffix", help="suffix to saved contigs that will be used to compute thresholds; empty by default", type=str)
    parser.set_defaults(etalon_contigs_prefix="")
    parser.set_defaults(etalon_contigs_suffix="")
    parser.add_argument("--force_to_update", help="comma-separated list of metrics for which threshold will be overwritten even if the metric satisfies current value; e.g. 'max_n50,min_n50,max_mis'")
    parser.set_defaults(force_to_update = "")
    parser.add_argument("--overwrite_if_satisfies", dest="overwrite_if_satisfies", help="overwrite threshold even if metric satisfies current value", action='store_true')
    parser.set_defaults(overwrite_if_satisfies=False)

    args = parser.parse_args()
    return args


BASIC_METRICS = ["N50", "LG50", "# misassemblies", "Genome fraction (%)",  "# indels per 100 kbp", "# mismatches per 100 kbp", "# local misassemblies"]


def metric_in_list(metric, mlist):
    return metric in mlist.split(',')


#use threshold update policy set by user
def shall_update_threshold(args, metric, entry):
    return (args.add_metrics == "all") or (args.add_metrics == "basic" and metric in BASIC_METRICS) or entry.assess or metric_in_list(entry.config_name, args.force_to_update)


#calculate new metric treshold
def calculate_treshold(entry, value):
    result = 0.0
    coeff = -1 if entry.should_be_higher_than_theshold else 1

    if entry.relative_delta:
        result = value * (1.0 + coeff * entry.delta_value)
    else:
        result = value + coeff * entry.delta_value

    if entry.is_int:
        result = math.floor(result) if entry.should_be_higher_than_theshold else math.ceil(result)
    result = max(0.0, result)
    return int(result) if entry.is_int else result


#update metrics using given quast report
def process_quast_report(args, report, limit_map):
    result_map = parse_report(report, limit_map)
    to_update = {}

    if args.overwrite_if_satisfies:
        log.log("Will overwrite thresholds even when metric satisfies")

    for metric in sorted(result_map.keys()):
        if metric in limit_map and len(limit_map[metric]) > 0:
            for entry in limit_map[metric]:
                param_name = entry.config_name

                if not shall_update_threshold(args, metric, entry):
                    log.log("Value for " + param_name + " will not be updated")
                    continue

                #that metric shouold be higher than threshold (e.g. N50)
                if entry.should_be_higher_than_theshold:
                    if result_map[metric] < entry.value or args.overwrite_if_satisfies or metric_in_list(param_name, args.force_to_update):
                        #either not satisies or overwrite anyways
                        to_update[param_name] = calculate_treshold(entry, result_map[metric])
                        log.log("New value for " + param_name + " is set to " + str(to_update[param_name]))
                    else:
                        log.log(param_name + " will not be updated, " + metric + " = " + str(result_map[metric]) + " >= " + str(entry.value))

                #that metric shouold be less than threshold (e.g. misassemblies)
                else:
                    if result_map[metric] > entry.value or args.overwrite_if_satisfies or metric_in_list(param_name, args.force_to_update):
                        #either not satisies or overwrite anyways
                        to_update[param_name] = calculate_treshold(entry, result_map[metric])
                        log.log("New value for " + param_name + " is set to " + str(to_update[param_name]))
                    else:
                        log.log(param_name + " will not be updated, " + metric + " = " + str(result_map[metric]) + " <= " + str(entry.value))
        else:
            log.log(metric + " = " + str(result_map[metric]))

    return to_update


# Run QUAST and assess its report for a single contig file
def quast_run_and_update(dataset_info, fn, output_dir, name, prefix, opts):
    if not os.path.exists(fn):
        log.err("File not found " + fn)
        return {}

    if 'quast_params' not in dataset_info.__dict__:
        log.log("quast_params are not set, will not run QUAST")
        return {}

    log.log("Processing " + fn)
    qcode = run_quast(dataset_info, [fn], output_dir, opts)
    if qcode != 0:
        log.err("Failed to run QUAST!")
        return {}

    limit_map = construct_rnaquast_limit_map(dataset_info, prefix, True) if dataset_info.mode == "rna" else construct_quast_limit_map(dataset_info, prefix, True)

    report_path = output_dir
    if dataset_info.mode == "meta":
        report_path = os.path.join(report_path, "combined_reference")
    if dataset_info.mode == "rna":
        report_path = os.path.join(report_path, "short_report.tsv")
    else:
        report_path = os.path.join(report_path, "report.tsv")

    return process_quast_report(args, report_path, limit_map)


def get_prefix(var):
    if var.startswith("sc_"):
        return "sc_"
    elif var.startswith("prelim_"):
        return "prelim_"
    return ""


#substitute and add params to config
def substitute_params(filename, var_dict):
    log.log(" === Substituting values in " + filename + " === ")
    lines = open(filename).readlines()
    vars_in_file = process_cfg.vars_from_lines(lines)

    for var, value in var_dict.items():
        var_prefix = get_prefix(var)
        assess_param = var_prefix + "assess"
        if assess_param not in vars_in_file:
            lines.extend(["\n", assess_param + " true", "\n"])
        vars_in_file = process_cfg.vars_from_lines(lines)

        if var not in vars_in_file:
            log.log(var + " is not in config, adding")
            assess_meta_info = vars_in_file[var_prefix + "assess"]
            lines = lines[:assess_meta_info.line_num + 1] + ["\n", assess_meta_info.indent + str(var) + " " + str(value) + "\n"] + lines[assess_meta_info.line_num + 1:]
            vars_in_file = process_cfg.vars_from_lines(lines)
        else:
            log.log("Substituting " + var)
            meta = vars_in_file[var]
            lines[meta.line_num] = meta.indent + str(var) + " " + str(value) + "\n"

    file = open(filename, "w")
    file.writelines(lines)
    file.close()


#run for given config
def update_thresholds_for_config(args, config):
    log.log("    >>>>>>>> Updating " + config)
    dataset_info = load_info(config)
    contigs = get_contigs_list(args, dataset_info)

    full_etalon_contigs_prefix = args.etalon_contigs_prefix
    etalon_contigs_suffix = args.etalon_contigs_suffix
    if  full_etalon_contigs_prefix == "":
        if 'contig_storage' in dataset_info.__dict__:
            full_etalon_contigs_prefix = os.path.join(dataset_info.contig_storage, "latest_")
            etalon_contigs_suffix = ""
        else:
            log.err("Etalon contigs are not set and contigs storage is not found in info file")
            return

    log.log("Using etalon contigs prefix: " + full_etalon_contigs_prefix)

    to_update = {}
    for name, file_name, prefix, opts, ext in contigs:
        if ext != "fasta" and ext != "fa":
            continue
        log.log("======= PROCESSING " + name.upper() + " =======")
        if prefix != "":
            prefix = prefix + "_"
        to_update.update(quast_run_and_update(dataset_info, full_etalon_contigs_prefix + name + etalon_contigs_suffix + "." + ext, os.path.join(args.tmp_dir, "QUAST_RESULTS_"  + dataset_info.name.upper() + "_" + name.upper()), name, prefix, opts))

    substitute_params(config, to_update)
    log.log("    <<<<<<<< Finished updating " + config)



#update threshold for all configs in folder
def update_thresholds_for_folder(args, config_dir):
    configs = sorted(glob.glob(os.path.join(config_dir, "*.info")))
    if len(configs) == 0:
        return

    for config in configs:
        update_thresholds_for_config(args, config)


#recursively update threshold for all configs in folder
def update_thresholds_recursive(args, base_dir):
    update_thresholds_for_folder(args, base_dir)
    
    files = sorted(glob.glob(os.path.join(base_dir, "*")))
    for f in files:
        if os.path.isdir(f):
            update_thresholds_recursive(args, f)


def check_args(args):
    correct = True

    if args.config == "" and args.config_dir == "":
        log.err("Either config file or config dir should be set")
        correct = False

    if args.config != "" and args.config_dir != "":
        log.err("Both config file or config dir cannot be set")
        correct = False

    if (args.etalon_contigs_prefix != "" or args.etalon_contigs_suffix) and args.config_dir != "":
        log.err("Etalon contigs cannot be set when using multiple config files from config")
        correct = False

    return correct


### main ###
try:
    if len(sys.argv) == 1:
        command = 'python {} -h'.format(sys.argv[0])
        subprocess.call(command, shell=True)
        sys.exit(1)

    sys.stderr = sys.stdout
    exit_code = 0
    args = parse_args_ut()
    if not check_args(args):
        log.err("Correct arguments and run again")
        sys.exit(3)

    if os.path.exists(args.tmp_dir) and not args.preserve_tmp_dir:
        shutil.rmtree(args.tmp_dir)
    if not os.path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir)

    if args.config != "":
        update_thresholds_for_config(args, args.config)
    else:
        update_thresholds_recursive(args, args.config_dir)

    if os.path.exists(args.tmp_dir) and not args.preserve_tmp_dir:
        shutil.rmtree(args.tmp_dir)

except SystemExit:
    raise

except:
    log.err("The following unexpected error occured during the run:")
    print_exc()
    sys.exit(239)
