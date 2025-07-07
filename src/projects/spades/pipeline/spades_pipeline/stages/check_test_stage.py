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

from ..options_storage import OptionStorage
options_storage = OptionStorage()
from . import stage
from ..commands_parser import Command


class CheckStageStage(stage.Stage):
    STAGE_NAME = "Check test"

    def get_command(self, cfg):
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "supplemetary", "check_test_script.py")]
        if options_storage.args.rna:
            args += ["--mode", "rna", "--result_transcripts_filename", self.output_files["result_transcripts_filename"]]
        else:
            if options_storage.args.plasmid:
                args += ["--mode", "plasmid"]
            else:
                args += ["--mode", "common"]

            args += ["--result_contigs_filename", self.output_files["result_contigs_filename"],
                     "--result_scaffolds_filename", self.output_files["result_scaffolds_filename"]]

        return [Command(STAGE=self.STAGE_NAME,
                                        path=sys.executable,
                                        args=args,
                                        short_name=self.short_name)]


def add_to_pipeline(pipeline, _, output_files, tmp_configs_dir, dataset_data, bin_home,
                    ext_python_modules_home, python_modules_home):
    if options_storage.args.test_mode:
        pipeline.add(CheckStageStage("check_test", output_files, tmp_configs_dir,
                                     dataset_data, bin_home, ext_python_modules_home, python_modules_home))
