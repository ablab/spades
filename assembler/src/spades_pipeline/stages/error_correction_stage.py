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

import commands_parser
import options_storage
from stages import stage
import process_cfg
import support
from process_cfg import merge_configs


class ECRunningToolStage(stage.Stage):
    def prepare_config_bh(self, filename, cfg, log):
        subst_dict = dict()
        subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset_yaml_filename)
        subst_dict["input_working_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
        subst_dict["output_dir"] = process_cfg.process_spaces(cfg.output_dir)
        subst_dict["general_max_iterations"] = options_storage.ITERATIONS
        subst_dict["general_max_nthreads"] = cfg.max_threads
        subst_dict["count_merge_nthreads"] = cfg.max_threads
        subst_dict["bayes_nthreads"] = cfg.max_threads
        subst_dict["expand_nthreads"] = cfg.max_threads
        subst_dict["correct_nthreads"] = cfg.max_threads
        subst_dict["general_hard_memory_limit"] = cfg.max_memory
        if "qvoffset" in cfg.__dict__:
            subst_dict["input_qvoffset"] = cfg.qvoffset
        if "count_filter_singletons" in cfg.__dict__:
            subst_dict["count_filter_singletons"] = cfg.count_filter_singletons
        if "read_buffer_size" in cfg.__dict__:
            subst_dict["count_split_buffer"] = cfg.read_buffer_size
        process_cfg.substitute_params(filename, subst_dict, log)

    def prepare_config_ih(self, filename, cfg, ext_python_modules_home):
        addsitedir(ext_python_modules_home)
        if sys.version.startswith("2."):
            import pyyaml2 as pyyaml
        elif sys.version.startswith("3."):
            import pyyaml3 as pyyaml
        data = pyyaml.load(open(filename))
        data["dataset"] = cfg.dataset_yaml_filename
        data["working_dir"] = cfg.tmp_dir
        data["output_dir"] = cfg.output_dir
        data["hard_memory_limit"] = cfg.max_memory
        data["max_nthreads"] = cfg.max_threads
        with open(filename, 'w') as f:
            pyyaml.dump(data, f,
                        default_flow_style=False, default_style='"', width=float("inf"))

    def generate_config(self, cfg):
        dst_configs = os.path.join(cfg.output_dir, "configs")
        if os.path.isdir(dst_configs):
            shutil.rmtree(dst_configs)
        if cfg.iontorrent:
            support.copy_tree(os.path.join(self.tmp_configs_dir, "ionhammer"), dst_configs, preserve_times=False)
            cfg_file_name = os.path.join(dst_configs, "ionhammer.cfg")
        else:
            support.copy_tree(os.path.join(self.tmp_configs_dir, "hammer"), dst_configs, preserve_times=False)
            cfg_file_name = os.path.join(dst_configs, "config.info")

        cfg.tmp_dir = support.get_tmp_dir(prefix="hammer_")
        if cfg.iontorrent:
            self.prepare_config_ih(cfg_file_name, cfg, self.ext_python_modules_home)
        else:
            self.prepare_config_bh(cfg_file_name, cfg, self.log)

    def get_command(self, cfg):
        dst_configs = os.path.join(cfg.output_dir, "configs")
        if cfg.iontorrent:
            cfg_file_name = os.path.join(dst_configs, "ionhammer.cfg")
        else:
            cfg_file_name = os.path.join(dst_configs, "config.info")

        if cfg.iontorrent:
            binary_name = "spades-ionhammer"
        else:
            binary_name = "spades-hammer"

        command = [commands_parser.Command(STAGE="Read error correction",
                                           path=os.path.join(self.bin_home, binary_name),
                                           args=[os.path.abspath(cfg_file_name)],
                                           config_dir=os.path.relpath(cfg.output_dir, options_storage.args.output_dir),
                                           short_name=self.short_name,
                                           del_after=[os.path.relpath(cfg.tmp_dir, options_storage.args.output_dir)],
                                           output_files=[self.output_files["corrected_dataset_yaml_filename"]])]
        return command


class ErrorCorrectionCompressingStage(stage.Stage):
    def get_command(self, cfg):
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "scripts", "compress_all.py"),
                "--input_file", self.output_files["corrected_dataset_yaml_filename"],
                "--ext_python_modules_home", self.ext_python_modules_home,
                "--max_threads", str(cfg.max_threads),
                "--output_dir", cfg.output_dir]
        if cfg.not_used_dataset_yaml_filename != "":
            args += ["--not_used_yaml_file", cfg.not_used_dataset_yaml_filename]
        if cfg.gzip_output:
            args.append("--gzip_output")

        command = [commands_parser.Command(STAGE="corrected reads compression",
                                           path=sys.executable,
                                           args=args,
                                           short_name=self.short_name)]
        return command


class ErrorCorrectionStage(stage.Stage):
    STAGE_NAME = "Read error correction"
    stages = []

    def __init__(self, cfg, *args):
        super(ErrorCorrectionStage, self).__init__(*args)

        self.cfg = merge_configs(cfg["error_correction"], cfg["common"])
        self.output_files["corrected_dataset_yaml_filename"] = os.path.join(self.cfg.output_dir, "corrected.yaml")

        self.cfg.not_used_dataset_yaml_filename = ""
        self.stages.append(ECRunningToolStage("ec_runtool",
                                              self.output_files, self.tmp_configs_dir, self.dataset_data,
                                              self.log, self.bin_home, self.ext_python_modules_home,
                                              self.python_modules_home))
        self.stages.append(
            ErrorCorrectionCompressingStage("ec_compress",
                                            self.output_files, self.tmp_configs_dir, self.dataset_data, self.log,
                                            self.bin_home, self.ext_python_modules_home,
                                            self.python_modules_home))


    def generate_config(self, cfg):
        if sys.version.startswith("2."):
            import pyyaml2 as pyyaml
        elif sys.version.startswith("3."):
            import pyyaml3 as pyyaml

        self.cfg = merge_configs(cfg["error_correction"], cfg["common"])
        self.output_files["corrected_dataset_yaml_filename"] = os.path.join(self.cfg.output_dir, "corrected.yaml")
        self.cfg.__dict__["dataset_yaml_filename"] = cfg["dataset"].yaml_filename

        addsitedir(self.ext_python_modules_home)

        if not os.path.isdir(self.cfg.output_dir):
            os.makedirs(self.cfg.output_dir)

        # not all reads need processing
        if support.get_lib_ids_by_type(self.dataset_data, options_storage.LONG_READS_TYPES):
            not_used_dataset_data = support.get_libs_by_type(self.dataset_data, options_storage.LONG_READS_TYPES)
            to_correct_dataset_data = support.rm_libs_by_type(self.dataset_data, options_storage.LONG_READS_TYPES)
            to_correct_dataset_yaml_filename = os.path.join(self.cfg.output_dir, "to_correct.yaml")
            self.cfg.not_used_dataset_yaml_filename = os.path.join(self.cfg.output_dir, "dont_correct.yaml")
            with open(to_correct_dataset_yaml_filename, 'w') as f:
                pyyaml.dump(to_correct_dataset_data, f,
                            default_flow_style=False, default_style='"', width=float("inf"))

            with open(self.cfg.not_used_dataset_yaml_filename, 'w') as f:
                pyyaml.dump(not_used_dataset_data, f,
                            default_flow_style=False, default_style='"', width=float("inf"))
            self.cfg.dataset_yaml_filename = to_correct_dataset_yaml_filename
        else:
            self.cfg.not_used_dataset_yaml_filename = ""

        for stage in self.stages:
            stage.generate_config(self.cfg)

    def get_command(self, cfg):
        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_start")] + \
               [x for stage in self.stages for x in stage.get_command(self.cfg)] + \
               [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_finish")]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log,
                    bin_home, ext_python_modules_home, python_modules_home):
    if "error_correction" in cfg:
        pipeline.add(ErrorCorrectionStage(cfg, "ec", output_files, tmp_configs_dir,
                                          dataset_data, log, bin_home, ext_python_modules_home, python_modules_home))
