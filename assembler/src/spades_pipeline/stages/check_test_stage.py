#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

import options_storage
from stages import stage
import commands_parser


class CheckStageStage(stage.Stage):
    STAGE_NAME = "Check test"

    def get_command(self, cfg):
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "scripts", "check_test_script.py")]
        if options_storage.args.truseq_mode:
            args += ["--mode", "truseq", "--truseq_long_reads_file", self.output_files["truseq_long_reads_file"]]
        elif options_storage.args.rna:
            args += ["--mode", "rna", "--result_transcripts_filename", self.output_files["result_transcripts_filename"]]
        else:
            if options_storage.args.plasmid:
                args += ["--mode", "plasmid"]
            else:
                args += ["--mode", "common"]

            args += ["--result_contigs_filename", self.output_files["result_contigs_filename"],
                     "--result_scaffolds_filename", self.output_files["result_scaffolds_filename"]]

        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path=sys.executable,
                                        args=args,
                                        short_name=self.short_name)]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log, bin_home,
                    ext_python_modules_home, python_modules_home):
    if options_storage.args.test_mode:
        pipeline.add(CheckStageStage("check_test", output_files, tmp_configs_dir,
                                     dataset_data, log, bin_home, ext_python_modules_home, python_modules_home))
