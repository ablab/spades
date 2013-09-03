#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import shutil
import support
import process_cfg
from site import addsitedir


def move_dataset_files(dataset_data, dst, ext_python_modules_home, max_threads, log, gzip=False):
    to_compress = []
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith('reads'):
                moved_reads_files = []
                for reads_file in value:
                    dst_filename = os.path.join(dst, os.path.basename(reads_file))
                    # TODO: fix problem with files with the same basenames in Hammer binary!
                    if not os.path.isfile(reads_file):
                        if (not gzip and os.path.isfile(dst_filename)) or (gzip and os.path.isfile(dst_filename + '.gz')):
                            support.warning('file with corrected reads (' + reads_file + ') is the same in several libraries', log)
                            if gzip:
                                dst_filename += '.gz'
                        else:
                            support.error('something went wrong and file with corrected reads (' + reads_file + ') is missing!', log)
                    else:
                        shutil.move(reads_file, dst_filename)
                        if gzip:
                            to_compress.append(dst_filename)
                            dst_filename += '.gz'
                    moved_reads_files.append(dst_filename)
                reads_library[key] = moved_reads_files
    if len(to_compress):
        pigz_path = support.which('pigz')
        if pigz_path:
            for reads_file in to_compress:
                support.sys_call([pigz_path, '-f', '-7', '-p', str(max_threads), reads_file], log)
        else:
            addsitedir(ext_python_modules_home)
            if sys.version.startswith('2.'):
                from joblib2 import Parallel, delayed
            elif sys.version.startswith('3.'):
                from joblib3 import Parallel, delayed
            n_jobs = min(len(to_compress), max_threads)
            outputs = Parallel(n_jobs=n_jobs)(delayed(support.sys_call)(['gzip', '-f', '-7', reads_file]) for reads_file in to_compress)
            for output in outputs:
                if output:
                    log.info(output)


def prepare_config_bh(filename, cfg, log):
    subst_dict = dict()

    subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset_yaml_filename)
    subst_dict["input_working_dir"] = process_cfg.process_spaces(os.path.abspath(cfg.tmp_dir))
    subst_dict["output_dir"] = process_cfg.process_spaces(os.path.abspath(cfg.tmp_dir))  # for only final corrected reads
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


def run_bh(result_filename, configs_dir, execution_home, cfg, ext_python_modules_home, log):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml

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

    command = [os.path.join(execution_home, "hammer"),
               os.path.abspath(cfg_file_name)]

    log.info("\n== Running read error correction tool: " + ' '.join(command) + "\n")
    support.sys_call(command, log)
    corrected_dataset_yaml_filename = os.path.join(cfg.tmp_dir, "corrected.yaml")
    if not os.path.isfile(corrected_dataset_yaml_filename):
        support.error("read error correction finished abnormally: " + corrected_dataset_yaml_filename + " not found!")
    corrected_dataset_data = pyyaml.load(open(corrected_dataset_yaml_filename, 'r'))
    if cfg.gzip_output:
        log.info("\n== Compressing corrected reads (with gzip)")
    move_dataset_files(corrected_dataset_data, cfg.output_dir, ext_python_modules_home, cfg.max_threads, log, cfg.gzip_output)
    corrected_dataset_yaml_filename = result_filename
    pyyaml.dump(corrected_dataset_data, open(corrected_dataset_yaml_filename, 'w'))
    log.info("\n== Dataset description file created: " + corrected_dataset_yaml_filename + "\n")

    shutil.rmtree(cfg.tmp_dir)
