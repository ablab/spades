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
from site import addsitedir

def prepare_config_bh(filename, cfg, log):
    subst_dict = dict()

    subst_dict["dataset"] = cfg.dataset_yaml_filename
    subst_dict["input_working_dir"] = os.path.abspath(cfg.tmp_dir)
    subst_dict["output_dir"] = os.path.abspath(cfg.tmp_dir)  # we will move only final corrected reads to output_dir
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


def run_bh(configs_dir, execution_home, cfg, ext_python_modules_home, log):
    addsitedir(ext_python_modules_home)
    import pyyaml

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
    corrected_dataset_yaml_filename = os.path.join(cfg.tmp_dir, "corrected.yaml")
    corrected_dataset_data = pyyaml.load(file(corrected_dataset_yaml_filename, 'r'))
    if cfg.gzip_output:
        log.info("\n== Compressing corrected reads (with gzip)")
    support.move_dataset_files(corrected_dataset_data, cfg.output_dir, log, cfg.gzip_output)
    corrected_dataset_yaml_filename = os.path.join(cfg.output_dir, "corrected.yaml")
    pyyaml.dump(corrected_dataset_data, file(corrected_dataset_yaml_filename, 'w'))
    log.info("\n== Dataset description file created: " + corrected_dataset_yaml_filename + "\n")

    shutil.rmtree(cfg.tmp_dir)
    return corrected_dataset_yaml_filename
