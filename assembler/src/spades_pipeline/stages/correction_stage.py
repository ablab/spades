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
from site import addsitedir

from stages import stage
import support
from process_cfg import merge_configs
import commands_parser
import options_storage

def prepare_config_corr(filename, cfg, ext_python_modules_home):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith("2."):
        import pyyaml2 as pyyaml
    elif sys.version.startswith("3."):
        import pyyaml3 as pyyaml
    data = pyyaml.load(open(filename))
    data["dataset"] = cfg.dataset
    data["output_dir"] = cfg.output_dir
    data["work_dir"] = cfg.tmp_dir
    # data["hard_memory_limit"] = cfg.max_memory
    data["max_nthreads"] = cfg.max_threads
    data["bwa"] = cfg.bwa
    with open(filename, 'w') as file_c:
        pyyaml.dump(data, file_c,
                    default_flow_style=False, default_style='"', width=float("inf"))


class CorrectionIterationStage(stage.Stage):
    def __init__(self, cfg, assembly_type, corrected, assembled, *args):
        super(CorrectionIterationStage, self).__init__(*args)
        self.assembly_type = assembly_type
        self.corrected = corrected
        self.assembled = assembled
        self.STAGE_NAME = "Mismatch correction %s" % assembly_type

        self.tmp_dir_for_corrector = os.path.join(cfg["common"].output_dir, "mismatch_corrector", self.assembly_type)
        cfg["mismatch_corrector"].__dict__["output_dir"] = self.tmp_dir_for_corrector
        self.cfg = merge_configs(cfg["mismatch_corrector"], cfg["common"])

    def get_command(self, cfg):
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "scripts", "correction_iteration_script.py"),
                "--corrected", self.corrected,
                "--assembled", self.assembled,
                "--assembly_type", self.assembly_type,
                "--output_dir", cfg["common"].output_dir,
                "--bin_home", self.bin_home]

        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path=sys.executable,
                                        args=args,
                                        config_dir=os.path.relpath(self.cfg.output_dir, options_storage.args.output_dir),
                                        short_name=self.short_name,
                                        del_after=[os.path.join(os.path.relpath(self.cfg.output_dir,
                                                                                options_storage.args.output_dir),
                                                                "tmp"),
                                                   os.path.relpath(self.cfg.tmp_dir, options_storage.args.output_dir)])]

    def generate_config(self, cfg):
        dst_configs = os.path.join(self.cfg.output_dir, "configs")
        if os.path.isdir(dst_configs):
            shutil.rmtree(dst_configs)
        support.copy_tree(os.path.join(self.tmp_configs_dir, "corrector"), dst_configs, preserve_times=False)
        cfg_file_name = os.path.join(dst_configs, "corrector.info")

        self.cfg.tmp_dir = support.get_tmp_dir(prefix="corrector_")
        prepare_config_corr(cfg_file_name, self.cfg, self.ext_python_modules_home)


class CorrectionStage(stage.Stage):
    stages = []
    STAGE_NAME = "Mismatch correction"

    def __init__(self, cfg, *args):
        super(CorrectionStage, self).__init__(*args)

        cfg["mismatch_corrector"].__dict__["dataset"] = cfg["dataset"].yaml_filename

        to_correct = dict()
        to_correct["contigs"] = \
            (self.output_files["result_contigs_filename"], self.output_files["assembled_contigs_filename"])
        to_correct["scaffolds"] = \
            (self.output_files["result_scaffolds_filename"], self.output_files["assembled_scaffolds_filename"])

        for assembly_type, (corrected, assembled) in to_correct.items():
            self.stages.append(CorrectionIterationStage(cfg, assembly_type, corrected, assembled,
                                                        "mc_%s" % assembly_type,
                                                        self.output_files,
                                                        self.tmp_configs_dir, self.dataset_data, self.log,
                                                        self.bin_home, self.ext_python_modules_home,
                                                        self.python_modules_home))

    def generate_config(self, cfg):
        for stage in self.stages:
            stage.generate_config(cfg)

    def get_command(self, cfg):
        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_start")] + \
               [x for stage in self.stages for x in stage.get_command(cfg)] + \
               [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_finish")]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log, bin_home,
                    ext_python_modules_home, python_modules_home):
    if "assembly" in cfg and "mismatch_corrector" in cfg:
        pipeline.add(CorrectionStage(cfg, "mc", output_files, tmp_configs_dir, dataset_data,
                                     log, bin_home, ext_python_modules_home, python_modules_home))
