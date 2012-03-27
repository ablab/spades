#!/usr/bin/env python

import os
import shutil
import sys
from os import path

spades_home = path.abspath(sys.path[0])

if os.path.realpath(__file__) == "/usr/bin/spades.py":
    spades_home = "/usr/share/spades/"

sys.path.append(os.path.join(spades_home, "src/tools/spades_pipeline/"))

import support
from process_cfg import *

def prepare_config(filename, cfg, prev_K, last_one):

    subst_dict = dict()
    cfg.output_dir = path.abspath(cfg.output_dir)
    cfg.input_dir = path.abspath(cfg.input_dir)

    subst_dict["dataset"]            = path.expandvars(cfg.dataset)
    subst_dict["input_dir"]          = path.dirname(subst_dict["dataset"])
    subst_dict["output_base"]        = cfg.build_path
    subst_dict["additional_contigs"] = path.join(cfg.build_path, "simplified_contigs.fasta")
    subst_dict["entry_point"]        = 'construction'

    if last_one:
        subst_dict["gap_closer_enable"] = "true"
        if cfg.paired_mode:
            subst_dict["paired_mode"] = "true"
    else:
        subst_dict["paired_mode"] = "false"
        subst_dict["gap_closer_enable"] = "false"

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
    os.environ["cfg"] = path.dirname(path.abspath(CONFIG_FILE))

    cfg = load_config_from_info_file(CONFIG_FILE)
    if cfg.has_key("bh"):   
        bh_cfg = cfg["bh"]
        if not bh_cfg.__dict__.has_key("output_dir"):
            bh_cfg.__dict__["output_dir"] = path.join(path.expandvars(cfg["common"].output_dir), "corrected")
        else:
            bh_cfg.__dict__["output_dir"] = path.expandvars(bh_cfg.output_dir)
        
        bh_cfg.__dict__["working_dir"] = path.join(bh_cfg.output_dir, "tmp")

        start_bh = True
        if path.exists(bh_cfg.output_dir):
            question = ["WARNING! Folder with corrected reads (" + bh_cfg.output_dir + ") already exists!", 
                        "Do you want to overwrite this folder and start error correction again? (y/n)"]
            answer = support.question_with_timer(question, 10)
            print 'answer is ', answer           
            if answer == 'n':
                start_bh = False
            else:
                shutil.rmtree(bh_cfg.output_dir)

        if start_bh:
            os.makedirs(bh_cfg.output_dir)
                
            log_filename = path.join(bh_cfg.output_dir, "bh.log")
            bh_cfg.__dict__["log_filename"] = log_filename

            print("\n== Error correction log can be found here: " + bh_cfg.log_filename + "\n")

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
                cfg["spades"].__dict__["dataset"] = run_bh(bh_cfg)
            except support.spades_error as err:
                print err.err_str
                err_code = err.code

            sys.stdout = old_stdout
            sys.stderr = old_stderr

            print("\n== Error correction finished. Log can be found here: " + bh_cfg.log_filename + "\n")
            if err_code:
               exit(err_code) 

    if cfg.has_key("spades"):
        spades_cfg = cfg["spades"]

        def build_folder(cfg):
            import datetime
            suffix = datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
            return path.join(path.expandvars(cfg.output_dir), 'spades_' + suffix)

        spades_cfg.__dict__["build_path"] = build_folder(cfg["common"])
        os.makedirs(spades_cfg.build_path)

        log_filename = path.join(spades_cfg.build_path, "spades.log")
        spades_cfg.__dict__["log_filename"] = log_filename

        print("\n== Assembling log can be found here: " + spades_cfg.log_filename + "\n")

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
        if err_code:
            exit(err_code)
    
    print("\n== SPAdes pipeline finished\n")

def prepare_config_bh(filename, cfg):

    subst_dict = dict()

    subst_dict["input_numfiles"]           = len(cfg.input_reads)
  
    for i in xrange(len(cfg.input_reads)):
        subst_dict["input_file_" + str(i)]  = path.expandvars(cfg.input_reads[i])
   
    subst_dict["input_gzipped"]             = bool_to_str(cfg.input_reads[0].endswith(".gz"))       
    subst_dict["input_working_dir"]         = cfg.working_dir
    subst_dict["general_max_iterations"]    = cfg.max_iterations
    subst_dict["general_max_nthreads"]      = cfg.max_threads  # TODO other nthreads values
    subst_dict["general_hard_memory_limit"] = cfg.max_memory
    # subst_dict["general_gzip"]              = bool_to_str(cfg.gzip_output)   # actually it is not for output (see hammer sources)

    substitute_params(filename, subst_dict)    

def run_bh(cfg):

    import build
    build.build_hammer(cfg, spades_home)
    
    dst_configs = path.join(cfg.output_dir, "configs")
    shutil.copytree(path.join(spades_home, "configs"), dst_configs)
    cfg_file_name = path.join(dst_configs, "hammer", "config.info")
        
    prepare_config_bh(cfg_file_name, cfg)

    command = path.join(os.getenv('HOME'), '.spades/precompiled/build_hammer', "hammer", "hammer") + " " + cfg_file_name

    support.sys_call(command)
    # support.sys_call("cp -r /home/alex/biolab/TEST/tmp " + cfg.working_dir)

    import bh_aux
    dataset_str = bh_aux.generate_dataset(cfg, cfg.working_dir)
    datasets_filename = path.abspath(path.join(cfg.output_dir, "dataset.info"))
    datasets_file = open(datasets_filename, "w")
    datasets_file.write(dataset_str)    
    datasets_file.close()
    
    shutil.rmtree(cfg.working_dir)

    return datasets_filename

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
