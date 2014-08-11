#!/usr/bin/python -O

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate coverage from raw file
import glob
import gzip
import sys
import os
import datetime
import getopt
from site import addsitedir
import logging
import shutil

from math import pow
from support import universal_sys_call, error
import options_storage

#profile = []
#insertions = {}
config = {}
#total_contigs
import os
import sys
import glob
import shutil
import support
import process_cfg
from site import addsitedir
from distutils import dir_util



def prepare_config_corr(filename, cfg, ext_python_modules_home):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml
    data = pyyaml.load(open(filename, 'r'))
    data["dataset"] = cfg.dataset
    data["output_dir"] = cfg.output_dir
    data["work_dir"] = cfg.output_dir + '/tmp'
    #data["hard_memory_limit"] = cfg.max_memory
    data["max_nthreads"] = cfg.max_threads
    data["bwa"] = cfg.bwa
    file_c = open(filename, 'w')
    pyyaml.dump(data, file_c)
    file_c.close()


def run_with_logger(to_correct, joblib_path, log=None, config_file=None):


    addsitedir(joblib_path)


    if not log:
        log = logging.getLogger('spades')
        log.setLevel(logging.DEBUG)

        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(logging.Formatter('%(message)s'))
        console.setLevel(logging.DEBUG)
        log.addHandler(console)

        log_filename = os.path.join(config["output_dirpath"], "corrector.log")
        log_handler = logging.FileHandler(log_filename, mode='w')
        log.addHandler(log_handler)
    path_to_bin = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../bin/corrector')
    path_to_config = os.path.join(os.path.dirname(os.path.realpath(__file__)) , '../../configs/corrector/corrector.info.template')
    if config_file:
        path_to_config = config_file
    run_str = path_to_bin + ' ' + path_to_config + ' ' + to_correct
    support.sys_call(run_str, log)



def run_corrector(corrected_dataset_yaml_filename, configs_dir, execution_home, cfg,
                ext_python_modules_home, log, to_correct):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml

    cfg.dataset_yaml_filename = corrected_dataset_yaml_filename
    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    dir_util.copy_tree(os.path.join(configs_dir, "corrector"), dst_configs, preserve_times=False)
    cfg_file_name = os.path.join(dst_configs, "corrector.info")
    # removing template configs
    for root, dirs, files in os.walk(dst_configs):
        for cfg_file in files:
            cfg_file = os.path.join(root, cfg_file)
            if cfg_file.endswith('.template'):
                if os.path.isfile(cfg_file.split('.template')[0]):
                    os.remove(cfg_file)
                else:
                    os.rename(cfg_file, cfg_file.split('.template')[0])

    cfg.tmp_dir = support.get_tmp_dir(prefix="corrector_")

    prepare_config_corr(cfg_file_name, cfg, ext_python_modules_home)
    binary_name = "corrector"

    command = [os.path.join(execution_home, binary_name),
               os.path.abspath(cfg_file_name)]

    log.info("\n== Running contig polishing tool: " + ' '.join(command) + "\n")


    log.info("\n== Dataset description file was created: " + cfg_file_name + "\n")
    run_with_logger(to_correct, ext_python_modules_home, log, cfg_file_name)




