#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

class Stage(object):
    def __init__(self, short_name, output_files, tmp_configs_dir,
                 dataset_data, bin_home, ext_python_modules_home, python_modules_home):
        self.short_name = short_name
        self.output_files = output_files
        self.tmp_configs_dir = tmp_configs_dir
        self.dataset_data = dataset_data
        self.bin_home = bin_home
        self.ext_python_modules_home = ext_python_modules_home
        self.python_modules_home = python_modules_home

    def generate_config(self, cfg):
        pass

    def get_command(self, _):
        return []
