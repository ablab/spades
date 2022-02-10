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

    # shutil.copyfile does not copy any metadata (time and permission), so one
    # cannot expect preserve_mode = False and preserve_times = True to work.
    def copy_tree(src, dst, preserve_times = True, preserve_mode = True):
        if preserve_mode == False:
            copy_fn = shutil.copyfile
        else:
            copy_fn = shutil.copy2

        # shutil.copytree preserves the timestamp, so we must update it afterwards.
        shutil.copytree(src, dst, copy_function = copy_fn, dirs_exist_ok = True)

        if preserve_times == False:
            for dirpath, _, filenames in os.walk(dst):
                os.utime(dirpath, None)
                for file in filenames:
                    os.utime(os.path.join(dirpath, file), None)
