#!/usr/bin/env python

import os
import shutil
import sys
from os import path

spades_home = "./"

if os.path.realpath(__file__) == "/usr/bin/spades.py":
    spades_home = "/usr/share/spades/"

sys.path.append(spades_home + "src/tools/spades_pipeline/")

import support
from process_cfg import *

def prepare_config(filename, cfg, prev_K, last_one):

    subst_dict = dict()
    cfg.output_dir = path.abspath(cfg.output_dir)
    cfg.input_dir = path.abspath(cfg.input_dir)

    subst_dict["dataset"]           = cfg.dataset
    subst_dict["input_dir"]         = path.dirname(cfg.dataset)
    subst_dict["output_base"]       = cfg.build_path
    subst_dict["additional_contigs"]= path.join(cfg.build_path, "simplified_contigs.fasta")
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

    CONFIG_FILE = spades_home + "spades_config.info"

    if os.path.isfile("spades_config.info") :
        CONFIG_FILE = "spades_config.info"

    if len(sys.argv) > 1 :
        if os.path.isfile(sys.argv[1]):
            CONFIG_FILE = sys.argv[1]
        else:
            print("Usage :")
            print("   ./spades.py <config file>")
            return

    print("Using config file: " + CONFIG_FILE)

    cfg = load_config_from_info_file(CONFIG_FILE)
    if cfg.has_key("bh"):    
        #TODO
        print "TODO: BayesHammer"

    if cfg.has_key("spades"):
#        spades_cfg = load_config_from_file(CONFIG_FILE) ##
        spades_cfg = cfg["spades"]

        def build_folder(cfg):
            import datetime
            suffix = datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
            return path.join(cfg.output_dir, 'spades_' + suffix)

        spades_cfg.__dict__["build_path"] = build_folder(cfg["common"])
        os.makedirs(spades_cfg.build_path)

        log_filename = path.join(spades_cfg.build_path, "spades.log")
        spades_cfg.__dict__["log_filename"] = log_filename

        print("\n== Log can be found here: " + spades_cfg.log_filename + "\n")

        log_file = open(log_filename, "w")

        old_stdout = sys.stdout
        old_stderr = sys.stderr

        if cfg["common"].output_to_console:
            sys.stderr = support.redirected_stream(log_file, sys.stderr)
            sys.stdout = support.redirected_stream(log_file, sys.stdout)
        else:
            sys.stderr = support.redirected_stream(log_file, None)
            sys.stdout = support.redirected_stream(log_file, None)

        err_code = 0
        try:
            run(spades_cfg)
        except support.spades_error as err:
            print err.err_str
            err_code = err.code

        sys.stdout = old_stdout
        sys.stderr = old_stderr

        print("\n== Assembling finished. Log can be found here: " + spades_cfg.log_filename + "\n")
        exit(err_code)
    
    print("\n== SPAdes pipeline finished\n")

def run(cfg):

    if type(cfg.iterative_K) is int:
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    import build
    build.build_spades_n_copy(cfg, spades_home)

    count = 0
    prev_K = None

    for K in cfg.iterative_K:
        count += 1

        dst_configs = path.join(cfg.build_path, str(K), "configs")
        shutil.copytree(path.join(spades_home, "configs"), dst_configs)
        cfg_file_name = path.join(dst_configs, "debruijn", "config.info")

        prepare_config(cfg_file_name, cfg, prev_K, count == len(cfg.iterative_K))
        prev_K = K

        execution_home = path.join(os.getenv('HOME'), '.spades/precompiled/build' + str(K))
        command = path.join(execution_home, "debruijn", "debruijn") + " " + path.abspath(cfg_file_name)

        print("\n== Running assembler: " + command + "\n")

        support.sys_call(command, execution_home)

        latest = path.join(cfg.build_path, "K%d" % (K), "latest")
        latest = os.readlink(latest)
        latest = os.path.join(cfg.build_path, "K%d" % (K), latest)
        os.symlink(os.path.relpath(latest, cfg.build_path), os.path.join(cfg.build_path, "link_K%d" % (K)))

    support.copy(os.path.join(latest, "result.info"), cfg.build_path)
    result = load_config_from_file(os.path.join(cfg.build_path, "result.info"))
    support.copy(result.contigs, cfg.build_path)

    print("\n== Running quality assessment tools: " + cfg.log_filename + "\n")
    cmd = "python " + path.join(spades_home, "src/tools/quality/quality.py") + " " + result.contigs
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
    print "All the resulting information can be found here: " + cfg.build_path
    print " * Resulting contigs are called " + os.path.split(result.contigs)[1]
    print " * Assessment of their quality is in " + qr + "/"
    print ""
    print "Thank you for using SPAdes!"

if __name__ == '__main__':
    main()
