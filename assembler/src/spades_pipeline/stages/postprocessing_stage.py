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


class PostprocessingStage(stage.Stage):
    STAGE_NAME = "Postprocessing"

    def get_command(self, cfg):
        # args: result_scaffolds_filename, assembled_scaffolds_filename, bin_home, ext_python_modules_home, output_dir, truseq_long_reads_file_base, dataset_yaml_file
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "scripts", "postprocessing_script.py"),
                "--result_scaffolds_filename", self.output_files["result_scaffolds_filename"],
                "--assembled_scaffolds_filename", self.output_files["assembled_scaffolds_filename"],
                "--bin_home", self.bin_home,
                "--ext_python_modules_home", self.ext_python_modules_home,
                "--output_dir", cfg["common"].output_dir,
                "--truseq_long_reads_file_base", self.output_files["truseq_long_reads_file_base"],
                "--dataset_yaml_file", options_storage.args.dataset_yaml_filename,
                "--threads", str(options_storage.args.threads)]

        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path=sys.executable,
                                        args=args,
                                        short_name=self.short_name)]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log,
                    bin_home, ext_python_modules_home, python_modules_home):
    if "assembly" in cfg and cfg["run_truseq_postprocessing"]:
        pipeline.add(
            PostprocessingStage("tpp", output_files, tmp_configs_dir, dataset_data, log,
                                bin_home, ext_python_modules_home, python_modules_home))
