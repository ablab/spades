#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import glob
import shutil
import support
import process_cfg
from site import addsitedir


def compress_dataset_files(dataset_data, ext_python_modules_home, max_threads, log):
    log.info("\n== Compressing corrected reads (with gzip)")
    to_compress = []
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith('reads'):
                compressed_reads_filenames = []
                for reads_file in value:
                    if not os.path.isfile(reads_file):
                        support.error('something went wrong and file with corrected reads (' + reads_file + ') is missing!', log)
                    to_compress.append(reads_file)
                    compressed_reads_filenames.append(reads_file + ".gz")
                reads_library[key] = compressed_reads_filenames
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


def remove_not_corrected_reads(output_dir):
    for not_corrected in glob.glob(os.path.join(output_dir, "*.bad.fastq")):
        os.remove(not_corrected)


def prepare_config_bh(filename, cfg, log):
    subst_dict = dict()

    subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset_yaml_filename)
    subst_dict["input_working_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
    subst_dict["output_dir"] = process_cfg.process_spaces(cfg.output_dir)
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


def prepare_config_ih(filename, cfg, ext_python_modules_home):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml

    data = pyyaml.load(open(filename, 'r'))
    data["dataset"] = process_cfg.process_spaces(cfg.dataset_yaml_filename)
    data["working_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
    data["output_dir"] = process_cfg.process_spaces(cfg.output_dir)
    data["hard_memory_limit"] = cfg.max_memory
    data["max_nthreads"] = cfg.max_threads
    pyyaml.dump(data, open(filename, 'w'))


def run_hammer(corrected_dataset_yaml_filename, configs_dir, execution_home, cfg,
               not_used_dataset_data, ext_python_modules_home, log):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml

    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    if cfg.iontorrent:
        shutil.copytree(os.path.join(configs_dir, "ionhammer"), dst_configs)
        cfg_file_name = os.path.join(dst_configs, "ionhammer.cfg")
    else:
        shutil.copytree(os.path.join(configs_dir, "hammer"), dst_configs)
        cfg_file_name = os.path.join(dst_configs, "config.info")
    # removing template configs
    for root, dirs, files in os.walk(dst_configs):
        for cfg_file in files:
            cfg_file = os.path.join(root, cfg_file)
            if cfg_file.endswith('.template'):
                if os.path.isfile(cfg_file.split('.template')[0]):
                    os.remove(cfg_file)
                else:
                    os.rename(cfg_file, cfg_file.split('.template')[0])

    cfg.tmp_dir = support.get_tmp_dir(prefix="hammer_")
    if cfg.iontorrent:
        prepare_config_ih(cfg_file_name, cfg, ext_python_modules_home)
        binary_name = "ionhammer"
    else:
        prepare_config_bh(cfg_file_name, cfg, log)
        binary_name = "hammer"

    command = [os.path.join(execution_home, binary_name),
               os.path.abspath(cfg_file_name)]

    log.info("\n== Running read error correction tool: " + ' '.join(command) + "\n")
    support.sys_call(command, log)
    if not os.path.isfile(corrected_dataset_yaml_filename):
        support.error("read error correction finished abnormally: " + corrected_dataset_yaml_filename + " not found!")
    corrected_dataset_data = pyyaml.load(open(corrected_dataset_yaml_filename, 'r'))
    remove_not_corrected_reads(cfg.output_dir)
    is_changed = False
    if cfg.gzip_output:
        is_changed = True
        compress_dataset_files(corrected_dataset_data, ext_python_modules_home, cfg.max_threads, log)
    if not_used_dataset_data:
        is_changed = True
        corrected_dataset_data += not_used_dataset_data
    if is_changed:
        pyyaml.dump(corrected_dataset_data, open(corrected_dataset_yaml_filename, 'w'))
    log.info("\n== Dataset description file was created: " + corrected_dataset_yaml_filename + "\n")

    if os.path.isdir(cfg.tmp_dir):
        shutil.rmtree(cfg.tmp_dir)