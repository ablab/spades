#!/usr/bin/env python

import os
import shutil
import sys

spades_home = os.path.abspath(sys.path[0])

if os.path.realpath(__file__) == "/usr/bin/spades.py":
    spades_home = "/usr/share/spades/"

sys.path.append(os.path.join(spades_home, "src/tools/spades_pipeline/"))
sys.path.append(os.path.join(spades_home, "src/tools/quality/"))
sys.path.append(os.path.join(spades_home, "src/tools/quality/libs"))

import support
from process_cfg import *

def error(err_str, prefix="== Error == "):
    print("\n\n" + prefix + " " + err_str + "\n\n")
    exit(1)

def warning(warn_str, prefix="== Warning == "):
    print("\n\n" + prefix + " " + warn_str + "\n\n")

def prepare_config_bh(filename, cfg):

    subst_dict = dict()
    cfg.working_dir = os.path.abspath(cfg.working_dir)
   
    #TODO single reads!

    subst_dict["input_numfiles"]           = len(cfg.paired_reads)

    for i, input_read in enumerate(cfg.paired_reads):
        subst_dict["input_file_" + str(i)]  = input_read

    subst_dict["input_gzipped"]             = bool_to_str(False)
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
    cfg.working_dir = os.path.abspath(cfg.working_dir)

    subst_dict["dataset"]            = cfg.dataset
    subst_dict["output_base"]        = cfg.working_dir
    subst_dict["additional_contigs"] = cfg.additional_contigs
    subst_dict["entry_point"]        = "construction"
    subst_dict["developer_mode"]     = str(cfg.developer_mode).lower()
    subst_dict["SAM_writer_enable"]  = bool_to_str(cfg.generate_sam_files)
    subst_dict["project_name"]       = ""

    if last_one:
        subst_dict["gap_closer_enable"] = "true"
        if cfg.paired_mode:
            subst_dict["paired_mode"] = "true"
    else:
        subst_dict["paired_mode"] = "false"
        subst_dict["gap_closer_enable"] = "false"

    if prev_K:
        subst_dict["use_additional_contigs"] = "true"
    else:
        subst_dict["use_additional_contigs"] = "false"

    substitute_params(filename, subst_dict)

def check_config(cfg, config_filename):

    if not cfg.has_key("dataset"):
        error("wrong config! You should specify 'dataset' section!")
        return False

    if (not cfg.has_key("error_correction")) and (not cfg.has_key("assembly")):
        error("wrong config! You should specify either 'error_correction' section (for reads error correction) or 'assembly' one (for assembling) or both!")
        return False

    if not cfg["common"].__dict__.has_key("output_dir"):
        error("wrong config! You should specify output_dir!")
        return False

    if not cfg["common"].__dict__.has_key("output_to_console"):
        cfg["common"].__dict__["output_to_console"] = True

    if not cfg["common"].__dict__.has_key("developer_mode"):
        cfg["common"].__dict__["developer_mode"] = False

    if not cfg["common"].__dict__.has_key("project_name"):
        cfg["common"].__dict__["project_name"] = os.path.splitext(os.path.basename(config_filename))[0]

    if not cfg["common"].__dict__.has_key("compilation_dir"):
        cfg["common"].__dict__["compilation_dir"] = os.path.join(os.getenv('HOME'), '.spades/precompiled/')
    else:
        cfg["common"].compilation_dir = os.path.abspath(cfg["common"].compilation_dir)

    cfg["common"].output_dir = os.path.join(os.path.abspath(os.path.expandvars(cfg["common"].output_dir)), cfg["common"].project_name)

    return True

def main():

    CONFIG_FILE = ""
    if os.path.isfile("spades_config.info"):
        CONFIG_FILE = "spades_config.info"
    elif os.path.isfile(os.path.join(spades_home, "spades_config.info")):
        CONFIG_FILE = os.path.join(spades_home, "spades_config.info")
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
       CONFIG_FILE = sys.argv[1]
    if not CONFIG_FILE:
        print >> sys.stderr, "Usage : python ", sys.argv[0], " <config file>"
        return

    print("Using config file: " + CONFIG_FILE)
    os.environ["cfg"] = os.path.dirname(os.path.abspath(CONFIG_FILE))

    cfg = load_config_from_info_file(CONFIG_FILE)

    if not check_config(cfg, CONFIG_FILE):
        return

    bh_dataset_filename = ""
    if cfg.has_key("error_correction"):
        bh_cfg = cfg["error_correction"]

        if not bh_cfg.__dict__.has_key("output_dir"):
            bh_cfg.__dict__["output_dir"] = os.path.join(os.path.expandvars(cfg["common"].output_dir), "corrected")
        else:
            bh_cfg.__dict__["output_dir"] = os.path.expandvars(bh_cfg.output_dir)
            warning("output_dir (" + bh_cfg.output_dir + ") will be used for error correction instead of the common one (" + cfg["common"].output_dir + ")")
        
        bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])
        
        bh_cfg.__dict__["working_dir"] = os.path.join(bh_cfg.output_dir, "tmp")

        bh_cfg.__dict__["dataset"] = os.path.join(bh_cfg.output_dir, cfg["common"].project_name + ".info")

        start_bh = True
        if os.path.exists(bh_cfg.output_dir):
            if os.path.exists(bh_cfg.dataset):
                question = ["WARNING! It looks like error correction was already done!", 
                            "Folder with corrected dataset " + bh_cfg.output_dir + " already exists!",
                            "Do you want to overwrite this folder and start error correction again?"]
                answer = support.question_with_timer(question, 10, 'n')
                if answer == 'n':
                    start_bh = False
                    bh_dataset_filename = bh_cfg.dataset
                    print("\n===== Error correction skipped\n")
                else:
                    os.remove(bh_cfg.dataset)
            #else:            
            #    shutil.rmtree(bh_cfg.output_dir)

        if start_bh:
            if not os.path.exists(bh_cfg.working_dir):
                os.makedirs(bh_cfg.working_dir)

            log_filename = os.path.join(bh_cfg.output_dir, "correction.log")
            bh_cfg.__dict__["log_filename"] = log_filename

            shutil.copy(CONFIG_FILE, bh_cfg.output_dir)

            # parsing dataset section
            bh_cfg.__dict__["single_cell"]  = cfg["dataset"].single_cell            
            bh_cfg.__dict__["paired_reads"] = []
            bh_cfg.__dict__["single_reads"] = []
            import bh_aux
            for key, value in cfg["dataset"].__dict__.iteritems():
                if type(value) != list:
                    value = [value]                
                if key.startswith("single_reads"):  
                    for item in value:
                        item = os.path.abspath(os.path.expandvars(item))
                        item = bh_aux.ungzip_if_needed(item, bh_cfg.working_dir)
                        bh_cfg.single_reads.append(item)
                elif key.startswith("paired_reads"):
                    cur_paired_reads = []
                    if len(value) == 1:
                        item = os.path.abspath(os.path.expandvars(value[0]))
                        cur_paired_reads = split_paired_file(item, output_folder)
                    elif len(value) == 2:
                        for item in value:
                            item = os.path.abspath(os.path.expandvars(item))
                            item = bh_aux.ungzip_if_needed(item, bh_cfg.working_dir)
                            cur_paired_reads.append(item)

                    if len(bh_cfg.paired_reads) == 0:
                        bh_cfg.paired_reads = cur_paired_reads
                    else:
                        bh_cfg.paired_reads = bh_aux.merge_paired_files(cur_paired_reads, bh_cfg.paired_reads, bh_cfg.working_dir)                          

            print("\n===== Error correction started. Log can be found here: " + bh_cfg.log_filename + "\n")

            tee = support.Tee(log_filename, 'w', console=bh_cfg.output_to_console)

            err_code = 0
            try:
                bh_dataset_filename = run_bh(bh_cfg)
            except support.spades_error, err:
                print err.err_str
                err_code = err.code

            tee.free()

            print("\n===== Error correction finished. Log can be found here: " + bh_cfg.log_filename + "\n")
            if err_code:
               exit(err_code)

    if cfg.has_key("assembly"):
        spades_cfg = merge_configs(cfg["assembly"], cfg["common"])        
        if not spades_cfg.__dict__.has_key("generate_sam_files"):
            spades_cfg.__dict__["generate_sam_files"] = False        

        def make_working_dir(output_dir):
            import datetime
            name = "spades_" + datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
            working_dir = os.path.join(output_dir, name)
            os.makedirs(working_dir)
            latest = os.path.join(output_dir, "latest")
            if os.path.islink(latest):
                os.remove(latest)
            if not os.path.exists(latest):
                os.symlink(name, latest)
            return working_dir

        spades_cfg.__dict__["working_dir"] = make_working_dir(spades_cfg.output_dir)
        spades_cfg.__dict__["log_filename"] = os.path.join(spades_cfg.working_dir, "assembly.log")
        spades_cfg.__dict__["result_contigs"] = os.path.join(spades_cfg.working_dir, spades_cfg.project_name + ".fasta")
        spades_cfg.__dict__["additional_contigs"] = os.path.join(spades_cfg.working_dir, "simplified_contigs.fasta")

        shutil.copy(CONFIG_FILE, spades_cfg.working_dir)      
        
        # dataset created during error correction
        if bh_dataset_filename:
            spades_cfg.__dict__["dataset"] = bh_dataset_filename
                
        if not spades_cfg.__dict__.has_key("dataset"):
            # creating dataset
            dataset_filename = os.path.join(spades_cfg.working_dir, cfg["common"].project_name + ".info")
            dataset_file = open(dataset_filename, 'w')
            for key, value in cfg["dataset"].__dict__.iteritems():
                dataset_file.write(key + '\t')

                if type(value) == bool:
                    dataset_file.write(bool_to_str(value))
                else:
                    dataset_file.write('"')
                    if type(value) != list:
                        value = [value]
                    for item in value:
                        item = os.path.abspath(os.path.expandvars(item))
                        dataset_file.write(str(item) + ' ')
                    dataset_file.write('"')

                dataset_file.write('\n')
            dataset_file.close()
            spades_cfg.__dict__["dataset"] = dataset_filename
        else:
            spades_cfg.dataset = os.path.abspath(os.path.expandvars(spades_cfg.dataset))
            shutil.copy(spades_cfg.dataset, spades_cfg.working_dir)

        print("\n===== Assembling started. Log can be found here: " + spades_cfg.log_filename + "\n")

        tee = support.Tee(spades_cfg.log_filename, 'w', console=spades_cfg.output_to_console)
        
        err_code = 0
        try:
            if cfg.has_key("quality_assessment"):
                run_spades(spades_cfg, cfg["quality_assessment"])
            else:
                run_spades(spades_cfg)
        except support.spades_error, err:
            print err.err_str
            err_code = err.code

        tee.free()

        print("\n===== Assembling finished. Log can be found here: " + spades_cfg.log_filename + "\n")
        if err_code:
            exit(err_code)

    print("\n===== SPAdes pipeline finished\n") 

def run_bh(cfg):

    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    shutil.copytree(os.path.join(spades_home, "configs"), dst_configs)
    cfg_file_name = os.path.join(dst_configs, "hammer", "config.info")

    prepare_config_bh(cfg_file_name, cfg)

    import build
    build.build_hammer(cfg, spades_home)

    execution_home = os.path.join(cfg.compilation_dir, 'build_hammer')
    command = os.path.join(execution_home, "hammer", "hammer") + " " + os.path.abspath(cfg_file_name)

    print("\n== Running error correction tool: " + command + "\n")
    support.sys_call(command)

    import bh_aux
    dataset_str = bh_aux.generate_dataset(cfg, cfg.working_dir, cfg.paired_reads)        
    dataset_filename = cfg.dataset
    dataset_file = open(dataset_filename, "w")
    dataset_file.write(dataset_str)
    dataset_file.close()
    print("\n== Dataset info-file created: " + dataset_filename + "\n")

    shutil.rmtree(cfg.working_dir)

    return dataset_filename

def run_spades(cfg, quality_cfg = None):

    if type(cfg.iterative_K) is int:
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    import build
    build.build_spades_n_copy(cfg, spades_home)

    count = 0
    prev_K = None

    for K in cfg.iterative_K:
        count += 1

        dst_configs = os.path.join(cfg.working_dir, "config_K" + str(K))
        os.mkdir(dst_configs)
        dst_configs = os.path.join(dst_configs, "configs")
        shutil.copytree(os.path.join(spades_home, "configs"), dst_configs)
        cfg_file_name = os.path.join(dst_configs, "debruijn", "config.info")

        prepare_config_spades(cfg_file_name, cfg, prev_K, count == len(cfg.iterative_K))
        prev_K = K

        execution_home = os.path.join(cfg.compilation_dir, 'build' + str(K))
        command = os.path.join(execution_home, "debruijn", "spades") + " " + os.path.abspath(cfg_file_name)

        print("\n== Running assembler: " + command + "\n")

        support.sys_call(command, execution_home)

        latest = os.path.join(cfg.working_dir, "K%d" % (K), "latest")
        latest = os.readlink(latest)
        latest = os.path.join(cfg.working_dir, "K%d" % (K), latest)
        #os.symlink(os.path.relpath(latest, cfg.working_dir), os.path.join(cfg.working_dir, "link_K%d" % (K)))
        os.symlink(latest, os.path.join(cfg.working_dir, "link_K%d" % (K)))  # python2.4 doesn't support os.path.relpath

    shutil.copyfile(os.path.join(latest, "final_contigs.fasta"), cfg.result_contigs)
    os.remove(cfg.additional_contigs)

    if quality_cfg:
        print("\n== Running quality assessment tools: " + cfg.log_filename + "\n")

        args = [cfg.result_contigs]
        
        if quality_cfg.__dict__.has_key("reference"):
            args.append("-R")
            args.append(os.path.abspath(os.path.expandvars(quality_cfg.reference)) )
        if quality_cfg.__dict__.has_key("genes"):
            args.append("-G")
            args.append(os.path.abspath(os.path.expandvars(quality_cfg.genes)) )
        if quality_cfg.__dict__.has_key("operons"):
            args.append("-O")
            args.append(os.path.abspath(os.path.expandvars(quality_cfg.operons)) )
        quality_output_dir = os.path.join(cfg.working_dir, "quality_results")
        args.append("-o")
        args.append(quality_output_dir)
        import quality
        quality.main(args, lib_dir=os.path.join(spades_home, "src/tools/quality/libs"))

    print ""
    print "All the resulting information can be found here: " + cfg.working_dir
    print " * Resulting contigs are called " + os.path.basename(cfg.result_contigs)
    if quality_cfg:
        print " * Assessment of their quality is in " + quality_output_dir + "/"
    print ""
    print "Thank you for using SPAdes!"

if __name__ == '__main__':
    main()
