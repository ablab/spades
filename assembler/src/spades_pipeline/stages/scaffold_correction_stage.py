#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
from distutils import dir_util

import commands_parser
import options_storage
from stages import stage
import process_cfg
from process_cfg import bool_to_str

# TODO copypast from iteraton stage
READS_TYPES_USED_IN_CONSTRUCTION = ["paired-end", "single", "hq-mate-pairs"]
READS_TYPES_USED_IN_RNA_SEQ = ["paired-end", "single", "trusted-contigs", "untrusted-contigs"]


def prepare_config_scaffold_correction(filename, cfg, log, saves_dir, K):
    subst_dict = dict()

    subst_dict["K"] = str(K)
    subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset)
    subst_dict["output_base"] = process_cfg.process_spaces(os.path.join(cfg.output_dir, "SCC"))
    subst_dict["tmp_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
    subst_dict["use_additional_contigs"] = bool_to_str(False)
    subst_dict["main_iteration"] = bool_to_str(False)
    subst_dict["entry_point"] = options_storage.BASE_STAGE
    subst_dict["load_from"] = saves_dir
    subst_dict["developer_mode"] = bool_to_str(cfg.developer_mode)
    subst_dict["max_threads"] = cfg.max_threads
    subst_dict["max_memory"] = cfg.max_memory

    # todo
    process_cfg.substitute_params(filename, subst_dict, log)


def add_configs(command, configs_dir):
    # Order matters here!
    mode_config_mapping = [("single_cell", "mda_mode"),
                           ("meta", "meta_mode"),
                           ("truseq_mode", "moleculo_mode"),
                           ("rna", "rna_mode"),
                           ("large_genome", "large_genome_mode"),
                           ("plasmid", "plasmid_mode")]
    # ("careful", "careful_mode"),
    for (mode, config) in mode_config_mapping:
        if options_storage.args.__dict__[mode]:
            if mode == "rna" or mode == "meta":
                command.append(os.path.join(configs_dir, "mda_mode.info"))
            command.append(os.path.join(configs_dir, config + ".info"))
    if options_storage.args.__dict__["careful"]:
        if options_storage.args.__dict__["single_cell"]:
            command.append(os.path.join(configs_dir, "careful_mda_mode.info"))
        else:
            command.append(os.path.join(configs_dir, "careful_mode.info"))

    # special case: extra config
    if options_storage.args.rna and options_storage.args.fast:
        command.append(os.path.join(configs_dir, "rna_fast_mode.info"))


class ScaffoldCorrectionStage(stage.Stage):
    def __init__(self, latest, *args):
        super(ScaffoldCorrectionStage, self).__init__(*args)
        self.latest = latest

    def generate_config(self, cfg):
        K = cfg.iterative_K[-1]
        latest = os.path.join(cfg.output_dir, "K%d" % K)
        K = options_storage.SCC_K
        data_dir = os.path.join(cfg.output_dir, "SCC", "K%d" % K)
        saves_dir = os.path.join(data_dir, "saves")
        dst_configs = os.path.join(data_dir, "configs")
        cfg_file_name = os.path.join(dst_configs, "config.info")

        if os.path.isdir(data_dir):
            shutil.rmtree(data_dir)
        os.makedirs(data_dir)

        dir_util.copy_tree(os.path.join(self.tmp_configs_dir, "debruijn"), dst_configs, preserve_times=False)

        scaffolds_file = os.path.join(latest, "scaffolds.fasta")
        if "read_buffer_size" in cfg.__dict__:
            construction_cfg_file_name = os.path.join(dst_configs, "construction.info")
            process_cfg.substitute_params(construction_cfg_file_name, {"read_buffer_size": cfg.read_buffer_size},
                                          self.log)
        process_cfg.substitute_params(os.path.join(dst_configs, "moleculo_mode.info"),
                                      {"scaffolds_file": scaffolds_file}, self.log)
        prepare_config_scaffold_correction(cfg_file_name, cfg, self.log, saves_dir, K)

    def get_command(self, cfg):
        K = options_storage.SCC_K
        data_dir = os.path.join(cfg.output_dir, "SCC", "K%d" % K)
        dst_configs = os.path.join(data_dir, "configs")
        cfg_file_name = os.path.join(dst_configs, "config.info")

        args = [cfg_file_name]
        add_configs(args, dst_configs)
        command = [commands_parser.Command(
            STAGE="SCC",
            path=os.path.join(self.bin_home, "spades-truseq-scfcorrection"),
            args=args,
            config_dir=os.path.relpath(data_dir, options_storage.args.output_dir),
            short_name=self.short_name)]
        return command
