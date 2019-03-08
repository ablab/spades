#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

class Stage(object):
    def __init__(self, short_name, output_files, tmp_configs_dir,
                 dataset_data, log, bin_home, ext_python_modules_home, python_modules_home):
        self.short_name = short_name
        self.output_files = output_files
        self.tmp_configs_dir = tmp_configs_dir
        self.dataset_data = dataset_data
        self.log = log
        self.bin_home = bin_home
        self.ext_python_modules_home = ext_python_modules_home
        self.python_modules_home = python_modules_home

    def generate_config(self, cfg):
        pass

    def get_command(self, cfg):
        return []
