#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

import options_storage
from stages import stage
import commands_parser


class TerminatingStage(stage.Stage):
    STAGE_NAME = "Terminate"

    def get_command(self, cfg):
        del_after = []
        if not cfg["common"].developer_mode:
            del_after.append(os.path.relpath(self.tmp_configs_dir, options_storage.args.output_dir))

        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name,
                                        del_after=del_after)]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log, bin_home,
                    ext_python_modules_home, python_modules_home):
    pipeline.add(TerminatingStage("terminate", output_files, tmp_configs_dir,
                                  dataset_data, log, bin_home,
                                  ext_python_modules_home, python_modules_home))