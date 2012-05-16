#!/usr/bin/env python

import os
import shutil
import sys
import spades_init
import getopt

spades_init.init()
spades_home = spades_init.spades_home

import support
from process_cfg import *

def error(err_str, prefix="== Error == "):
    print >>sys.stderr, "\n\n" + prefix + " " + err_str + "\n\n"
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

    subst_dict["dataset"]                = cfg.dataset
    subst_dict["output_base"]            = cfg.working_dir
    subst_dict["additional_contigs"]     = cfg.additional_contigs
    subst_dict["entry_point"]            = "construction"
    subst_dict["developer_mode"]         = bool_to_str(cfg.developer_mode)
    subst_dict["SAM_writer_enable"]      = bool_to_str(cfg.generate_sam_files and last_one)
    subst_dict["align_original_reads"]   = bool_to_str(cfg.align_original_reads)
    subst_dict["project_name"]           = ""
    subst_dict["gap_closer_enable"]      = bool_to_str(last_one and cfg.gap_closer)
    subst_dict["paired_mode"]            = bool_to_str(last_one and cfg.paired_mode)
    subst_dict["use_additional_contigs"] = bool_to_str(prev_K)

    substitute_params(filename, subst_dict)

def check_config(cfg, default_project_name=""):

    ## checking mandatory sections

    if not cfg.has_key("dataset") and not (cfg.has_key("assembly") and cfg["assembly"].__dict__.has_key("dataset")):
        error("wrong config! You should specify 'dataset' section!")
        return False

    if (not cfg.has_key("error_correction")) and (not cfg.has_key("assembly")):
        error("wrong config! You should specify either 'error_correction' section (for reads error correction) or 'assembly' one (for assembling) or both!")
        return False

    ## checking existence of all files in dataset section

    no_files_with_reads = True
    for k, v in cfg["dataset"].__dict__.iteritems():
        if k.find("reads") != -1:
            no_files_with_reads = False
            if type(v) != list:
                v = [v]
            for reads_file in v:
                if not os.path.isfile(os.path.expandvars(reads_file)):
                    error("file with reads doesn't exist! " + os.path.expandvars(reads_file))
                    return False
    
    if no_files_with_reads:
        error("wrong config! You should specify at least one file with reads!")
        return False
   
    ## setting default values if needed
   
    # common 
    if not cfg["common"].__dict__.has_key("output_dir"):
        cfg["common"].__dict__["output_dir"]     = 'spades_output'

    if not cfg["common"].__dict__.has_key("max_threads"):
        cfg["common"].__dict__["max_threads"]    = 16

    if not cfg["common"].__dict__.has_key("max_memory"):
        cfg["common"].__dict__["max_memory"]     = 250

    if not cfg["common"].__dict__.has_key("output_to_console"):
        cfg["common"].__dict__["output_to_console"] = True

    if not cfg["common"].__dict__.has_key("developer_mode"):
        cfg["common"].__dict__["developer_mode"]    = False

    if not cfg["common"].__dict__.has_key("project_name"):
        cfg["common"].__dict__["project_name"]      = default_project_name

    cfg["common"].__dict__["compilation_dir"] = os.path.join(os.getenv('HOME'), '.spades/precompiled/')    

    cfg["common"].output_dir = os.path.join(os.path.abspath(os.path.expandvars(cfg["common"].output_dir)), cfg["common"].project_name)

    # dataset
    if cfg.has_key("dataset"):
        if not cfg["dataset"].__dict__.has_key("single_cell"):
            cfg["dataset"].__dict__["single_cell"] = False            

    # error_correction
    if cfg.has_key("error_correction"):
        if not cfg["error_correction"].__dict__.has_key("max_iterations"):
            cfg["error_correction"].__dict__["max_iterations"] = 2    
        if not cfg["error_correction"].__dict__.has_key("gzip_output"):
            cfg["error_correction"].__dict__["gzip_output"] = True    

    # assembly
    if cfg.has_key("assembly"):
        if not cfg["assembly"].__dict__.has_key("iterative_K"):
            cfg["assembly"].__dict__["iterative_K"] = [21, 33, 55]
        if not cfg["assembly"].__dict__.has_key("generate_sam_files"):
            cfg["assembly"].__dict__["generate_sam_files"] = False        
        if not cfg["assembly"].__dict__.has_key("gap_closer"):
            cfg["assembly"].__dict__["gap_closer"] = True
            
    return True

long_options = "project-name= paired= paired1= paired2= single= k-mers= sc help".split()
short_options = "n:p:1:2:s:k:Sh"

def check_file(f, message=''):
    if not os.path.isfile(f):
        error("file not found (%s): %s" % (message, f))
    return f

def usage():
    print >>sys.stderr, 'SPAdes genome assembler'
    print >>sys.stderr, 'Usage:'
    print >>sys.stderr, sys.argv[0], ' <config file>'
    print >>sys.stderr, 'or'
    print >>sys.stderr, sys.argv[0], ' [options described below] -n <project name>'
    print >>sys.stderr, ""
    print >>sys.stderr, "Options with parameters:"
    print >>sys.stderr, "-n\t--project-name\tName of the project [Mandatory parameter]"
    print >>sys.stderr, "-p\t--paired\tFile with interlaced left and right paired end reads"
    print >>sys.stderr, "-1\t--paired1\tFile with left paired end reads"
    print >>sys.stderr, "-2\t--paired2\tFile with right paired end reads"
    print >>sys.stderr, "-s\t--single\tFile with unpaired reads"
    print >>sys.stderr, "-k\t--k-mers\tComma-separated list of odd values for k-mer (vertex) sizes. Default is 21,33,55"    
    print >>sys.stderr, ""
    print >>sys.stderr, "Options without parameters:"
    print >>sys.stderr, "-S\t--sc\tShould be set if input data was obtained with MDA (single-cell) technology"
    print >>sys.stderr, "-h\t--help\tPrint this usage message"    

def main():

    CONFIG_FILE = ""
    options = None

    if len(sys.argv) == 1:
        if os.path.isfile("spades_config.info"):
            CONFIG_FILE = "spades_config.info"
        elif os.path.isfile(os.path.join(spades_home, "spades_config.info")):
            CONFIG_FILE = os.path.join(spades_home, "spades_config.info")
    elif len(sys.argv) == 2:
        if os.path.isfile(sys.argv[1]):
            CONFIG_FILE = sys.argv[1]   
    else:
        try:
            options, not_options = getopt.gnu_getopt(sys.argv, short_options, long_options)
        except getopt.GetoptError, err:
            print >>sys.stderr, err
            print >>sys.stderr
            usage()
            sys.exit(1)

    if not CONFIG_FILE and not options:
        usage()
        sys.exit(1)

    # all parameters are stored here
    cfg = dict()

    if CONFIG_FILE:
        cfg = load_config_from_info_file(CONFIG_FILE)

        if not check_config(cfg, os.path.splitext(os.path.basename(CONFIG_FILE))[0] ):
            return
    
    elif options:
        project_name = ""
        paired  = []
        paired1 = []
        paired2 = []
        single  = []
        k_mers  = []
        single_cell = False

        for opt, arg in options:        
            if opt in ('-n', "--project-name"):
                project_name = arg
            elif opt in ('-p', "--paired"):
                paired.append(check_file(arg, 'paired'))
            elif opt in ('-1', "--paired1"):
                paired1.append(check_file(arg, 'left paired'))
            elif opt in ('-2', "--paired2"):
                paired2.append(check_file(arg, 'right paired'))
            elif opt in ('-s', "--single"):
                single.append(check_file(arg, 'single'))
            elif opt in ('-k', "--k-mers"):
                k_mers = map(int, arg.split(","))
            elif opt in ('-S', "--sc"):
                single_cell = True
            elif opt in ('-h', "--help"):
                usage()
                sys.exit(0)
            else:
                raise ValueError    

        if not project_name:
            error("the project name is not set! It is a mandatory parameter.")

        if len(paired1) != len(paired2):
            error("the number of files with left paired reads is not equal to the number of files with right paired reads!")           

        if not paired and not paired1 and not single:
            error("you should specify either paired reads or single reads or both!")
        
        # filling cfg
        cfg["common"]           = load_config_from_vars(dict())
        cfg["dataset"]          = load_config_from_vars(dict()) 
        cfg["error_correction"] = load_config_from_vars(dict())      
        cfg["assembly"]         = load_config_from_vars(dict())

        # filling reads
        paired_counter = 0
        if paired:
            for read in paired:
                paired_counter += 1
                cfg["dataset"].__dict__["paired_reads" + str(paired_counter)] = read
    
        if paired1:
            for i in range(len(paired1)):
                paired_counter += 1
                cfg["dataset"].__dict__["paired_reads" + str(paired_counter)] = [paired1[i], paired2[i]]

        if single:
            cfg["dataset"].__dict__["single_reads"] = single

        # filling other parameters
        if k_mers:
            cfg["assembly"].__dict__["iterative_K"] = k_mers

        if single_cell:
            cfg["dataset"].__dict__["single_cell"] = True

        if not check_config(cfg, project_name):
            return


    print("\n======= SPAdes pipeline started\n")

    if CONFIG_FILE:
        print("Using config file: " + CONFIG_FILE)
        os.environ["cfg"] = os.path.dirname(os.path.abspath(CONFIG_FILE))   

    bh_dataset_filename = ""
    if cfg.has_key("error_correction"):
        bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])
      
        bh_cfg.output_dir = os.path.join(os.path.expandvars(bh_cfg.output_dir), "corrected")
        
        if bh_cfg.__dict__.has_key("tmp_dir"):
            bh_cfg.__dict__["working_dir"] = os.path.expandvars(bh_cfg.tmp_dir)
        else:
            bh_cfg.__dict__["working_dir"] = os.path.join(bh_cfg.output_dir, "tmp")

        bh_cfg.__dict__["dataset"] = os.path.join(bh_cfg.output_dir, cfg["common"].project_name + ".dataset")

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
                    shutil.rmtree(bh_cfg.output_dir)

        if start_bh:            

            if not os.path.exists(bh_cfg.working_dir):
                os.makedirs(bh_cfg.working_dir)

            log_filename = os.path.join(bh_cfg.output_dir, "correction.log")
            bh_cfg.__dict__["log_filename"] = log_filename

            print("\n===== Error correction started. Log can be found here: " + bh_cfg.log_filename + "\n")
            tee = support.Tee(log_filename, 'w', console=bh_cfg.output_to_console)

            if CONFIG_FILE:
                shutil.copy(CONFIG_FILE, bh_cfg.output_dir)
            
            # parsing dataset section
            bh_cfg.__dict__["single_cell"]  = cfg["dataset"].single_cell            
            bh_cfg.__dict__["paired_reads"] = []
            bh_cfg.__dict__["single_reads"] = []
            import bh_aux
            for key, value in cfg["dataset"].__dict__.iteritems():
                if type(value) != list:
                    value = [value]                

                # saving original reads to dataset
                if key.find("_reads") != -1:
                    quoted_value = '"'
                    for item in value:
                        quoted_value += os.path.abspath(os.path.expandvars(item)) + ' '
                    quoted_value += '"'
                    bh_cfg.__dict__["original_" + key] = quoted_value
            
                # saving reference to dataset in developer_mode
                if bh_cfg.developer_mode:
                    if cfg.has_key("quality_assessment") and cfg["quality_assessment"].__dict__.has_key("reference"):
                        bh_cfg.__dict__["reference_genome"] = os.path.abspath(os.path.expandvars(cfg["quality_assessment"].reference)) 

                if key.startswith("single_reads"):  
                    for item in value:
                        item = os.path.abspath(os.path.expandvars(item))
                        item = bh_aux.ungzip_if_needed(item, bh_cfg.working_dir)
                        if len(bh_cfg.single_reads) == 0:
                            bh_cfg.single_reads.append(item)
                        else:
                            bh_cfg.single_reads[0] = bh_aux.merge_single_files(item, bh_cfg.single_reads[0], bh_cfg.working_dir)

                elif key.startswith("paired_reads"):
                    cur_paired_reads = []
                    if len(value) == 1:
                        item = os.path.abspath(os.path.expandvars(value[0]))
                        cur_paired_reads = bh_aux.split_paired_file(item, bh_cfg.working_dir)
                    elif len(value) == 2:
                        for item in value:
                            item = os.path.abspath(os.path.expandvars(item))
                            item = bh_aux.ungzip_if_needed(item, bh_cfg.working_dir)
                            cur_paired_reads.append(item)

                    if len(bh_cfg.paired_reads) == 0:
                        bh_cfg.paired_reads = cur_paired_reads
                    else:
                        bh_cfg.paired_reads = bh_aux.merge_paired_files(cur_paired_reads, bh_cfg.paired_reads, bh_cfg.working_dir)                                      
            
            bh_dataset_filename = run_bh(bh_cfg)
            
            tee.free()
            print("\n===== Error correction finished. Log can be found here: " + bh_cfg.log_filename + "\n")            

    result_contigs_filename = ""
    if cfg.has_key("assembly"):        
        spades_cfg = merge_configs(cfg["assembly"], cfg["common"])        
        
        if cfg.has_key("error_correction"):
            spades_cfg.__dict__["align_original_reads"] = True
        else:
            spades_cfg.__dict__["align_original_reads"] = False
        if cfg["dataset"].__dict__.has_key("paired_reads") or cfg["dataset"].__dict__.has_key("paired_reads.1"):
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
        spades_cfg.__dict__["log_filename"] = os.path.join(spades_cfg.working_dir, "assembly.log")
        spades_cfg.__dict__["result_contigs"] = os.path.join(spades_cfg.working_dir, spades_cfg.project_name + ".fasta")
        spades_cfg.__dict__["additional_contigs"] = os.path.join(spades_cfg.working_dir, "simplified_contigs.fasta")
        final_contigs_folder = os.path.join(spades_cfg.output_dir, "contigs")

        make_link(os.path.basename(spades_cfg.working_dir), os.path.join(spades_cfg.output_dir, "latest"))
        if os.path.exists(final_contigs_folder):
            shutil.rmtree(final_contigs_folder)

        print("\n===== Assembling started. Log can be found here: " + spades_cfg.log_filename + "\n")
        tee = support.Tee(spades_cfg.log_filename, 'w', console=spades_cfg.output_to_console)

        if CONFIG_FILE:
            shutil.copy(CONFIG_FILE, spades_cfg.working_dir)      
        
        # dataset created during error correction
        if bh_dataset_filename:
            spades_cfg.__dict__["dataset"] = bh_dataset_filename
                
        if not spades_cfg.__dict__.has_key("dataset"):
            # creating dataset
            dataset_filename = os.path.join(spades_cfg.working_dir, cfg["common"].project_name + ".dataset")
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
            # saving reference to dataset in developer_mode
            if spades_cfg.developer_mode:
                if cfg.has_key("quality_assessment") and cfg["quality_assessment"].__dict__.has_key("reference"):
                    dataset_file.write("reference_genome" + '\t')
                    dataset_file.write(os.path.abspath(os.path.expandvars(cfg["quality_assessment"].reference)) + '\n')

            dataset_file.close()
            spades_cfg.__dict__["dataset"] = dataset_filename
        else:
            spades_cfg.dataset = os.path.abspath(os.path.expandvars(spades_cfg.dataset))
            shutil.copy(spades_cfg.dataset, spades_cfg.working_dir)      
            # cheat for developers: if dataset was set in 'assembly' section we should use reference from it
            if cfg.has_key("quality_assessment") and not bh_dataset_filename:
                for key_to_del in ["reference", "genes", "operons"]:
                    if cfg["quality_assessment"].__dict__.has_key(key_to_del):
                        del cfg["quality_assessment"].__dict__[key_to_del]

                dataset_cfg = load_config_from_file(spades_cfg.dataset)
                for k,v in dataset_cfg.__dict__.iteritems():
                    if k == "reference_genome":
                        cfg["quality_assessment"].__dict__["reference"] = os.path.join(os.path.dirname(spades_cfg.dataset), v)
                
        
        
        result_contigs_filename = run_spades(spades_cfg)              
        
        tee.free()
        print("\n===== Assembling finished. Log can be found here: " + spades_cfg.log_filename + "\n")        

        make_link(os.path.basename(spades_cfg.working_dir), os.path.join(spades_cfg.output_dir, "latest_success"))
        os.makedirs(final_contigs_folder)
        shutil.copy(result_contigs_filename, final_contigs_folder)
        shutil.copy(spades_cfg.log_filename, final_contigs_folder)

    quality_final_report = ""
    if cfg.has_key("quality_assessment") and result_contigs_filename:
        
        quality_cfg = merge_configs(cfg["quality_assessment"], cfg["common"])

        quality_cfg.__dict__["working_dir"] = os.path.dirname(result_contigs_filename)
        quality_cfg.__dict__["log_filename"] = os.path.join(quality_cfg.working_dir, "quality.log")
        quality_cfg.__dict__["result_contigs"] = result_contigs_filename        

        print("\n===== Quality assessment started. Log can be found here: " + quality_cfg.log_filename + "\n")        
        tee = support.Tee(quality_cfg.log_filename, 'w', console=quality_cfg.output_to_console)
        
        quality_final_report = run_quality(quality_cfg)            
        
        tee.free()
        print("\n===== Quality assessment finished. Log can be found here: " + quality_cfg.log_filename + "\n")        

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

    import build
    build.build_hammer(cfg, spades_home)

    execution_home = os.path.join(cfg.compilation_dir, 'build_hammer')
    command = os.path.join(execution_home, "hammer", "hammer") + " " + os.path.abspath(cfg_file_name)

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
    if cfg.developer_mode:
        before_RR_contigs = os.path.join(os.path.dirname(cfg.result_contigs), "contigs_before_RR.fasta")
        shutil.copyfile(os.path.join(latest, "contigs_before_RR.fasta"), before_RR_contigs)
    os.remove(cfg.additional_contigs) 

    return cfg.result_contigs

def run_quality(cfg):        

    args = [cfg.result_contigs]

    if cfg.developer_mode:
        before_RR_contigs = os.path.join(os.path.dirname(cfg.result_contigs), "contigs_before_RR.fasta")
        args.append(before_RR_contigs)
    
    if cfg.__dict__.has_key("reference"):
        args.append("-R")
        args.append(os.path.abspath(os.path.expandvars(cfg.reference)) )
    if cfg.__dict__.has_key("genes"):
        args.append("-G")
        args.append(os.path.abspath(os.path.expandvars(cfg.genes)) )
    if cfg.__dict__.has_key("operons"):
        args.append("-O")
        args.append(os.path.abspath(os.path.expandvars(cfg.operons)) )
    quality_output_dir = os.path.join(cfg.working_dir, "quality_results")
    args.append("-o")
    args.append(quality_output_dir)
    import quality
    quality.main(args, lib_dir=os.path.join(spades_home, "src/tools/quality/libs"))
    
    return os.path.join(quality_output_dir, "quality.txt")

if __name__ == '__main__':
    main()
