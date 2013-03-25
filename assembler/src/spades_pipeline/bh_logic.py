#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import shutil

import support
import process_cfg

def prepare_config_bh(filename, cfg, log):
    subst_dict = dict()

    subst_dict["dataset"] = cfg.dataset_yaml_filename
    subst_dict["input_working_dir"] = os.path.abspath(cfg.tmp_dir)
    subst_dict["output_dir"] = os.path.abspath(cfg.output_dir)
    subst_dict["general_max_iterations"] = cfg.max_iterations
    subst_dict["general_max_nthreads"] = cfg.max_threads
    subst_dict["count_merge_nthreads"] = cfg.max_threads
    subst_dict["bayes_nthreads"] = cfg.max_threads
    subst_dict["expand_nthreads"] = cfg.max_threads
    subst_dict["correct_nthreads"] = cfg.max_threads
    subst_dict["general_hard_memory_limit"] = cfg.max_memory

    if "qvoffset" in cfg.__dict__:
        subst_dict["input_qvoffset"] = cfg.qvoffset

    process_cfg.substitute_params(filename, subst_dict, log)


def run_bh(configs_dir, execution_home, cfg, log):
    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    shutil.copytree(os.path.join(configs_dir, "hammer"), dst_configs)
    cfg_file_name = os.path.join(dst_configs, "config.info")
    # removing template configs
    for root, dirs, files in os.walk(dst_configs):
        for cfg_file in files:
            cfg_file = os.path.join(root, cfg_file)
            if cfg_file.endswith('.info.template'):
                if os.path.isfile(cfg_file.split('.template')[0]):
                    os.remove(cfg_file)
                else:
                    os.rename(cfg_file, cfg_file.split('.template')[0])

    prepare_config_bh(cfg_file_name, cfg, log)

    command = os.path.join(execution_home, "hammer") + " " +\
               os.path.abspath(cfg_file_name)

    log.info("\n== Running read error correction tool: " + command + "\n")
    support.sys_call(command, log)

    import bh_aux
    dataset_str = bh_aux.generate_dataset(cfg)
    dataset_filename = cfg.dataset
    dataset_file = open(dataset_filename, "w")
    dataset_file.write(dataset_str)
    dataset_file.close()
    log.info("\n== Dataset description file created: " + dataset_filename + "\n")

    #TODO: uncomment removing when BH will generate corrected reads directly to <output>/corrected (NOT IN ../tmp)
    #shutil.rmtree(cfg.tmp_dir)

    return dataset_filename
