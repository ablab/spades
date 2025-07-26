#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

import pyyaml3 as pyyaml

from ..commands_parser import Command
from . import stage


#delete tmp files in ouput folder
class BeforeStartStage(stage.Stage):
    STAGE_NAME = "Before start"
    stages = []

    def __init__(self, cfg, *args):
        super(BeforeStartStage, self).__init__(*args)
        output_dir = cfg["common"].output_dir
        self.tmp_files = []

        if (os.path.isfile(os.path.join(output_dir, "run_spades.yaml"))):
            previous_pipeline = pyyaml.load(open(os.path.join(output_dir, "run_spades.yaml")))
            for previous_stage in previous_pipeline:
                self.tmp_files += previous_stage["del_after"]

    def get_command(self, cfg):
        return [Command(stage=self.STAGE_NAME,
                        path="true",
                        args=[],
                        short_name=self.short_name,
                        del_after=self.tmp_files)]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data,
                    bin_home, ext_python_modules_home, python_modules_home):
    pipeline.add(BeforeStartStage(cfg, "before_start", output_files, tmp_configs_dir,
                                  dataset_data, bin_home, ext_python_modules_home, python_modules_home))