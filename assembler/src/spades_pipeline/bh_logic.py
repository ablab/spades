#!/usr/bin/env python

import os
import shutil

import support
import process_cfg

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

    process_cfg.substitute_params(filename, subst_dict)


def run_bh(spades_home, execution_home, cfg):
    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    shutil.copytree(os.path.join(spades_home, "configs", "hammer"), dst_configs)
    cfg_file_name = os.path.join(dst_configs, "config.info")
    # removing template configs
    for root, dirs, files in os.walk(dst_configs):
        for cfg_file in files:
            if cfg_file.endswith('.template'):
                os.remove(os.path.join(root, cfg_file))

    prepare_config_bh(cfg_file_name, cfg)

    command = ""
    if "use_jemalloc" in cfg.__dict__ and os.path.isfile("jemalloc.sh"):
        command = os.path.abspath("jemalloc.sh") + " "

    command += os.path.join(execution_home, "hammer") + " " +\
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

    shutil.rmtree(cfg.tmp_dir)

    return dataset_filename