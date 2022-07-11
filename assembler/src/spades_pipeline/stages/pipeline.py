#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil

import support
import options_storage

class Pipeline(object):
    stages = []

    # copying configs before all computations (to prevent its changing at run time)
    def copy_configs(self, cfg, spades_home, tmp_configs_dir):
        if os.path.isdir(tmp_configs_dir):
            shutil.rmtree(tmp_configs_dir)
        if not os.path.isdir(tmp_configs_dir):
            if options_storage.args.configs_dir:
                support.copy_tree(options_storage.args.configs_dir, tmp_configs_dir, preserve_times=False,
                                   preserve_mode=False)
            else:
                support.copy_tree(os.path.join(spades_home, "configs"), tmp_configs_dir, preserve_times=False,
                                   preserve_mode=False)

    def add(self, stage):
        self.stages.append(stage)

    def get_commands(self, cfg):
        commands = []
        for stage in self.stages:
            commands += stage.get_command(cfg)
        return commands

    def generate_configs(self, cfg, spades_home, tmp_configs_dir):
        self.copy_configs(cfg, spades_home, tmp_configs_dir)
        for stage in self.stages:
            stage.generate_config(cfg)
