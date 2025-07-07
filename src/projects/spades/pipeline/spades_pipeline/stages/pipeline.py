#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil

from ..options_storage import OptionStorage
options_storage = OptionStorage()
from ..support import copy_tree


class Pipeline(object):
    stages = []

    # copying configs before all computations (to prevent its changing at run time)
    def copy_configs(self, config_dirs, tmp_configs_dir):
        if os.path.isdir(tmp_configs_dir):
            shutil.rmtree(tmp_configs_dir)
        if not os.path.isdir(tmp_configs_dir):
            if options_storage.args.configs_dir:
                copy_tree(options_storage.args.configs_dir, tmp_configs_dir, preserve_times=False,
                                   preserve_mode=False)
            else:
                for src_dir, target_dir in config_dirs:
                    copy_tree(src_dir, os.path.join(tmp_configs_dir, target_dir),
                              preserve_times=False, preserve_mode=False)

    def add(self, stage):
        self.stages.append(stage)

    def get_commands(self, cfg):
        commands = []
        for stage in self.stages:
            commands += stage.get_command(cfg)
        return commands

    def generate_configs(self, cfg, config_dirs, tmp_configs_dir):
        self.copy_configs(config_dirs, tmp_configs_dir)
        for stage in self.stages:
            stage.generate_config(cfg)
