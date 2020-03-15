#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

if sys.version.startswith("2."):
    import pyyaml2 as pyyaml
elif sys.version.startswith("3."):
    import pyyaml3 as pyyaml

import commands_parser
from stages import stage

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
        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name,
                                        del_after=self.tmp_files)]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log,
                    bin_home, ext_python_modules_home, python_modules_home):
    pipeline.add(BeforeStartStage(cfg, "before_start", output_files, tmp_configs_dir,
                                  dataset_data, log, bin_home, ext_python_modules_home, python_modules_home))