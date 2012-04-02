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

def error(err_str, prefix="== Error == "):
    print("\n\n" + prefix + " " + err_str + "\n\n")
    exit(1)

def warning(warn_str, prefix="== Warning == "):
    print("\n\n" + prefix + " " + warn_str + "\n\n")

def prepare_config_bh(filename, cfg):

    subst_dict = dict()
    cfg.working_dir = path.abspath(cfg.working_dir)

    import glob
    input_reads = []
    if not type(cfg.input_reads) is list:
        cfg.input_reads = [cfg.input_reads]
    for read in cfg.input_reads:
        input_reads.extend(glob.glob(path.abspath(path.expandvars(read))))
    cfg.input_reads = input_reads

    if len(cfg.input_reads) == 1:
        import bh_aux
        cfg.input_reads = bh_aux.split_paired_file(cfg)

    subst_dict["input_numfiles"]           = len(cfg.input_reads)
  
    for i in xrange(len(cfg.input_reads)):
        subst_dict["input_file_" + str(i)]  = cfg.input_reads[i]
   
    subst_dict["input_gzipped"]             = bool_to_str(cfg.input_reads[0].endswith(".gz"))       
    subst_dict["input_working_dir"]         = cfg.working_dir
    subst_dict["general_max_iterations"]    = cfg.max_iterations
    subst_dict["general_max_nthreads"]      = cfg.max_threads 
    subst_dict["count_merge_nthreads"]      = cfg.max_threads 
    subst_dict["bayes_nthreads"]            = cfg.max_threads 
    subst_dict["expand_nthreads"]           = cfg.max_threads 
    subst_dict["correct_nthreads"]          = cfg.max_threads                
    subst_dict["general_hard_memory_limit"] = cfg.max_memory
    
    if cfg.__dict__.has_key("qvoffset"):
        subst_dict["input_qvoffset"]        = cfg.qvoffset

    substitute_params(filename, subst_dict) 

def prepare_config_spades(filename, cfg, prev_K, last_one):

    subst_dict = dict()
    cfg.working_dir = path.abspath(cfg.working_dir)

    subst_dict["dataset"]            = cfg.dataset    
    subst_dict["output_base"]        = cfg.working_dir
    subst_dict["additional_contigs"] = path.join(cfg.working_dir, "simplified_contigs.fasta")
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

    substitute_params(filename, subst_dict)

def check_config(cfg, config_filename):

    if (not cfg.has_key("bh")) and (not cfg.has_key("spades")):
        error("wrong config! You should specify either 'bh' section (for reads error correction) or 'spades' one (for assembling) or both!")
        return False
    
    if not cfg["common"].__dict__.has_key("output_dir"):
        error("wrong config! You should specify output_dir!")
        return False

    if not cfg["common"].__dict__.has_key("output_to_console"):
        cfg["common"].__dict__["output_to_console"] = True        

    if not cfg["common"].__dict__.has_key("project_name"):
        cfg["common"].__dict__["project_name"] = path.splitext(path.basename(config_filename))[0]

    cfg["common"].output_dir = path.join(path.abspath(path.expandvars(cfg["common"].output_dir)), cfg["common"].project_name)

    return True

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

    if not check_config(cfg, CONFIG_FILE):
        return     

    created_dataset_filename = ""
    if cfg.has_key("bh"):   
        bh_cfg = cfg["bh"]
        if not bh_cfg.__dict__.has_key("output_dir"):
            bh_cfg.__dict__["output_dir"] = path.join(cfg["common"].output_dir, "corrected")
        else:
            bh_cfg.__dict__["output_dir"] = path.expandvars(bh_cfg.output_dir)
            warning("output_dir (" + bh_cfg.output_dir + ") will be used for error correction instead of the common one (" + cfg["common"].output_dir + ")")
        
        bh_cfg.__dict__["working_dir"] = path.join(bh_cfg.output_dir, "tmp")

        start_bh = True
        if path.exists(bh_cfg.output_dir):
            question = ["WARNING! Folder with corrected reads (" + bh_cfg.output_dir + ") already exists!", 
                        "Do you want to overwrite this folder and start error correction again?"]
            answer = support.question_with_timer(question, 10, 'n')
            if answer == 'n':
                start_bh = False
                print("\n===== Error correction skipped\n")   
            else:
                shutil.rmtree(bh_cfg.output_dir)

        if start_bh:
            os.makedirs(bh_cfg.working_dir)
                
            log_filename = path.join(bh_cfg.output_dir, "bh.log")
            bh_cfg.__dict__["log_filename"] = log_filename

            shutil.copy(CONFIG_FILE, bh_cfg.output_dir)

            print("\n===== Error correction started. Log can be found here: " + bh_cfg.log_filename + "\n")

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
                created_dataset_filename = run_bh(bh_cfg)                
            except support.spades_error as err:
                print err.err_str
                err_code = err.code

            sys.stdout = old_stdout
            sys.stderr = old_stderr

            print("\n===== Error correction finished. Log can be found here: " + bh_cfg.log_filename + "\n")
            if err_code:
               exit(err_code) 

    if cfg.has_key("spades"):
        spades_cfg = cfg["spades"]
        if not spades_cfg.__dict__.has_key("measure_quality"):
            spades_cfg.__dict__["measure_quality"] = False

        if created_dataset_filename != "":
            if cfg["spades"].__dict__.has_key("dataset"):
                warning("dataset created during error correction (" + created_dataset_filename + ") will be used instead of the one from config file (" + cfg["spades"].dataset + ")!")            
            cfg["spades"].__dict__["dataset"] = created_dataset_filename 
        spades_cfg.dataset = path.abspath(path.expandvars(spades_cfg.dataset)) 
        if not path.isfile(spades_cfg.dataset):
            error("dataset " + spades_cfg.dataset + " doesn't exist!")
            exit(1)
        if not check_dataset(spades_cfg.dataset):
            error("" "incorrect dataset " + spades_cfg.dataset + ": list of files with reads should be arounded with double quotes!")
            exit(1)

        def make_working_dir(cfg):
            import datetime
            name = "spades_" + datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
            output_dir  = cfg.output_dir
            working_dir = path.join(output_dir, name)
            os.makedirs(working_dir)
            latest = path.join(output_dir, "the_latest")
            if os.path.islink(latest):
                os.remove(latest)
            if not os.path.exists(latest):
                os.symlink(name, latest)
            return working_dir

        spades_cfg.__dict__["working_dir"] = make_working_dir(cfg["common"])

        log_filename = path.join(spades_cfg.working_dir, "spades.log")
        spades_cfg.__dict__["log_filename"] = log_filename

        shutil.copy(CONFIG_FILE, spades_cfg.working_dir)
        shutil.copy(spades_cfg.dataset, spades_cfg.working_dir)
        #spades_cfg.dataset = path.join(spades_cfg.working_dir, path.basename(spades_cfg.dataset))
        #correct_dataset(spades_cfg.dataset)

        print("\n===== Assembling started. Log can be found here: " + spades_cfg.log_filename + "\n")

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
            run_spades(spades_cfg)
        except support.spades_error as err:
            print err.err_str
            err_code = err.code

        sys.stdout = old_stdout
        sys.stderr = old_stderr

        print("\n===== Assembling finished. Log can be found here: " + spades_cfg.log_filename + "\n")
        if err_code:
            exit(err_code)
    
    print("\n===== SPAdes pipeline finished\n") 

def run_bh(cfg):    
    
    dst_configs = path.join(cfg.output_dir, "configs")
    shutil.copytree(path.join(spades_home, "configs"), dst_configs)
    cfg_file_name = path.join(dst_configs, "hammer", "config.info")
        
    prepare_config_bh(cfg_file_name, cfg)

    import build
    build.build_hammer(cfg, spades_home)

    execution_home = path.join(os.getenv('HOME'), '.spades/precompiled/build_hammer')
    command = path.join(execution_home, "hammer", "hammer") + " " + path.abspath(cfg_file_name)
    
    print("\n== Running BayesHammer: " + command + "\n")
    support.sys_call(command)

    import bh_aux
    dataset_str = bh_aux.generate_dataset(cfg, cfg.working_dir, cfg.input_reads)
    if not cfg.__dict__.has_key("dataset_name"):
        cfg.__dict__["dataset_name"] = "dataset"
    dataset_filename = path.abspath(path.join(cfg.output_dir, cfg.dataset_name + ".info"))
    dataset_file = open(dataset_filename, "w")
    dataset_file.write(dataset_str)    
    dataset_file.close()
    print("\n== Dataset created: " + dataset_filename + "\n")    

    shutil.rmtree(cfg.working_dir)

    return dataset_filename

def run_spades(cfg):

    if type(cfg.iterative_K) is int:
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    import build
    build.build_spades_n_copy(cfg, spades_home)

    count = 0
    prev_K = None

    for K in cfg.iterative_K:
        count += 1

        dst_configs = path.join(cfg.working_dir, "config_K" + str(K), "configs")
        shutil.copytree(path.join(spades_home, "configs"), dst_configs)
        cfg_file_name = path.join(dst_configs, "debruijn", "config.info")

        prepare_config_spades(cfg_file_name, cfg, prev_K, count == len(cfg.iterative_K))
        prev_K = K

        execution_home = path.join(os.getenv('HOME'), '.spades/precompiled/build' + str(K))
        command = path.join(execution_home, "debruijn", "debruijn") + " " + path.abspath(cfg_file_name)

        print("\n== Running assembler: " + command + "\n")

        support.sys_call(command, execution_home)

        #dataset_id = path.splitext(path.basename(cfg.dataset))[0]
        latest = path.join(cfg.working_dir, "K%d" % (K), "latest")
        latest = os.readlink(latest)
        latest = path.join(cfg.working_dir, "K%d" % (K), latest)
        os.symlink(os.path.relpath(latest, cfg.working_dir), os.path.join(cfg.working_dir, "link_K%d" % (K)))

    support.copy(os.path.join(latest, "final_contigs.fasta"), cfg.working_dir)
    result_contigs = os.path.join(cfg.working_dir, "final_contigs.fasta")

    if cfg.measure_quality:
        print("\n== Running quality assessment tools: " + cfg.log_filename + "\n")
        cmd = "python " + path.join(spades_home, "src/tools/quality/quality.py") + " " + result_contigs
        dataset_filename = cfg.dataset
        dataset = load_config_from_info_file(dataset_filename)["common"]

        if dataset.__dict__.has_key("reference_genome"):
            cmd += " -R " + path.join(path.dirname(dataset_filename), dataset.reference_genome)
        if dataset.__dict__.has_key("genes"):
            cmd += " -G " + path.join(path.dirname(dataset_filename), dataset.genes)
        if dataset.__dict__.has_key("operons"):
            cmd += " -O " + path.join(path.dirname(dataset_filename), dataset.operons)
        qr = "quality_results"
        cmd += " -o " + os.path.join(cfg.working_dir, qr)
        support.sys_call(cmd)

    print ""
    print "All the resulting information can be found here: " + cfg.working_dir
    print " * Resulting contigs are called " + os.path.split(result_contigs)[1]
    if cfg.measure_quality:
        print " * Assessment of their quality is in " + qr + "/"
    print ""
    print "Thank you for using SPAdes!"

if __name__ == '__main__':
    main()
