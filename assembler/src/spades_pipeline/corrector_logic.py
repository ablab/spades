#!/usr/bin/python -O

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import shutil
import support
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
    data["work_dir"] = os.path.join(cfg.output_dir, 'tmp')
    #data["hard_memory_limit"] = cfg.max_memory
    data["max_nthreads"] = cfg.max_threads
    data["bwa"] = cfg.bwa
    file_c = open(filename, 'w')
    pyyaml.dump(data, file_c)
    file_c.close()



def run_corrector(configs_dir, execution_home, cfg,
                ext_python_modules_home, log, to_correct, result):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml

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
               os.path.abspath(cfg_file_name), os.path.abspath(to_correct)]

    log.info("\n== Running contig polishing tool: " + ' '.join(command) + "\n")


    log.info("\n== Dataset description file was created: " + cfg_file_name + "\n")

    support.sys_call(command, log)
    if not os.path.isfile(result):
        support.error("Mismatch correction finished abnormally: " + result + " not found!")




