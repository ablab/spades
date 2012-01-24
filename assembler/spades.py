#!/usr/bin/env python

import os
import sys

sys.path.append("src/tools/spades_pipeline/") # is it OK to have relative path here? what if user runs it from other directory?

import support
from process_cfg import *

def prepare_config(filename, build_path, cfg, prev_K, last_one):

    subst_dict = dict()

    subst_dict["dataset"]           = cfg.dataset
    subst_dict["input_dir"]         = cfg.input_dir
    subst_dict["output_base"]       = cfg.output_dir
    subst_dict["additional_contigs"]= build_path + "simplified_contigs.fasta"

    if cfg.paired_mode and last_one:
        subst_dict["paired_mode"] = "true"
    else:
        subst_dict["paired_mode"] = "false"

    if prev_K is not None:
        subst_dict["use_additional_contigs"] = "true"
    else:
        subst_dict["use_additional_contigs"] = "false"

#    if cfg.measure_quality:
#        print("Quality measuring implementation is in progress")

    substitute_params(filename, subst_dict)

def main():
    if len(sys.argv) == 2:
        config_file_name = sys.argv[1]
    elif len(sys.argv) == 1:
        config_file_name = "spades_config.info"
    else:
        print("Usage ./spades.py <config file>")

    cfg = load_config_from_file(config_file_name)

    def build_folder(cfg):
        import datetime
        suffix = datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
        return cfg.output_dir + cfg.dataset + r"/build_" + suffix + r"/"

    build_path = build_folder(cfg)
    os.makedirs(build_path)

    log_filename = build_path + "spades.log"
    cfg.__dict__["log_filename"] = log_filename

    print("\n== Log can be found here: " + cfg.log_filename + "\n")

    log_file = open(log_filename, "w")

    old_stdout = sys.stdout
    old_stderr = sys.stderr

    if cfg.output_to_console:
        sys.stderr = support.redirected_stream(log_file, sys.stderr)
        sys.stdout = support.redirected_stream(log_file, sys.stdout)
    else:
        sys.stderr = support.redirected_stream(log_file, None)
        sys.stdout = support.redirected_stream(log_file, None)

    # --

    if type(cfg.iterative_K) is int:
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    import build
    build.build_spades_n_copy(cfg, build_path)

    count = 0
    prev_K = None

    for K in cfg.iterative_K:
        count += 1

        path = build_path + "/" + str(K) + "/"
        cfg_file_name = path + "configs/debruijn/config.info"
        prepare_config(cfg_file_name, build_path, cfg, prev_K, count == len(cfg.iterative_K))
        prev_K = K

        support.sys_call(cfg.output_to_console, cfg.log_filename, path + "debruijn " + cfg_file_name)

    sys.stdout = old_stdout
    sys.stderr = old_stderr

    print("\n== Assembling finished. Log can be found here: " + cfg.log_filename + "\n")

if __name__ == '__main__':
  main()
