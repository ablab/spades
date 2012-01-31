#!/usr/bin/env python

import os
import sys

sys.path.append("src/tools/spades_pipeline/") # is it OK to have relative path here? what if user runs it from other directory?

import support
from process_cfg import *

def prepare_config(filename, cfg, prev_K, last_one):

    subst_dict = dict()

    subst_dict["dataset"]           = cfg.dataset
    subst_dict["input_dir"]         = cfg.input_dir
    subst_dict["output_base"]       = cfg.output_dir
    subst_dict["additional_contigs"]= cfg.build_path + "simplified_contigs.fasta"
    subst_dict["entry_point"]       = 'construction'


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

    cfg.__dict__["build_path"] = build_folder(cfg)
    os.makedirs(cfg.build_path)

    log_filename = cfg.build_path + "spades.log"
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

    err_code = 0
    try:
        run(cfg)
    except support.spades_error as err:
        print err.err_str
        err_code = err.code

    sys.stdout = old_stdout
    sys.stderr = old_stderr

    print("\n== Assembling finished. Log can be found here: " + cfg.log_filename + "\n")
    exit(err_code)

def run(cfg):


    if type(cfg.iterative_K) is int:
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    import build
    build.build_spades_n_copy(cfg)

    count = 0
    prev_K = None

    for K in cfg.iterative_K:
        count += 1

        path = cfg.build_path + "/" + str(K) + "/"
        cfg_file_name = path + "configs/debruijn/config.info"
        prepare_config(cfg_file_name, cfg, prev_K, count == len(cfg.iterative_K))
        prev_K = K

        support.sys_call(os.path.join(path, "debruijn") + " " + cfg_file_name)
        latest = os.path.join(cfg.output_dir, cfg.dataset, "K%d" % (K), "latest")
        latest = os.readlink(latest)
        latest = os.path.join(cfg.output_dir, cfg.dataset, "K%d" % (K), latest)
        os.symlink(os.path.relpath(latest, build_path), os.path.join(build_path, "link_K%d" % (K)))


    fn = os.path.join(cfg.output_dir, cfg.dataset + "/K" + str(prev_K) + "/latest/result.info")
    result = load_config_from_file(fn)

    support.sys_call("cp " + result.contigs + " " + cfg.build_path)

    print("\n== Running quality assessment tools: " + cfg.log_filename + "\n")
    cmd = "python src/tools/quality/quality.py " + result.contigs
    if result.reference:
        cmd += " -R " + result.reference
#    if result.genes:
#        cmd += " -G " + result.genes
#    if result.operons:
#        cmd += " -O " + result.operons
    qr = "quality_results"
    cmd += " -o " + os.path.join(cfg.build_path, qr)
    support.sys_call(cmd)

    print ""
    print "All the resulting information can be found here: " + build_path
    print " * Resulting contigs are called " + os.path.split(result.contigs)[1]
    print " * Assessment of their quality is in " + qr + "/"
    print ""
    print "Thank you for using SPAdes!"

if __name__ == '__main__':
    main()
