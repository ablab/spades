#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
import sys

import options_storage
from stages import stage
import commands_parser


class BreakingScaffoldsStage(stage.Stage):
    STAGE_NAME = "Breaking scaffolds"

    def get_command(self, cfg):
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "scripts", "breaking_scaffolds_script.py"),
                "--result_scaffolds_filename", self.output_files["result_scaffolds_filename"],
                "--misc_dir", self.output_files["misc_dir"],
                "--threshold_for_breaking_scaffolds", str(options_storage.THRESHOLD_FOR_BREAKING_SCAFFOLDS)]

        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path=sys.executable,
                                        args=args,
                                        short_name=self.short_name)]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log,
                    bin_home, ext_python_modules_home, python_modules_home):
    pipeline.add(BreakingScaffoldsStage("bs", output_files,
                                        tmp_configs_dir, dataset_data, log, bin_home,
                                        ext_python_modules_home, python_modules_home))
