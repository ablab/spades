#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil

import commands_parser
import options_storage
from stages import stage
import support
import process_cfg
from process_cfg import bool_to_str


# FIXME double with scaffold correction stage
def add_configs(command, configs_dir, cfg):
    # Order matters here!
    mode_config_mapping = [("isolate", "isolate_mode"),
                           ("single_cell", "mda_mode"),
                           ("meta", "meta_mode"),
                           ("truseq_mode", "moleculo_mode"),
                           ("rna", "rna_mode"),
                           ("large_genome", "large_genome_mode"),
                           ("plasmid", "plasmid_mode"),
                           ("metaviral", "metaviral_mode"),
                           ("metaplasmid", "metaplasmid_mode"),
                           ("rnaviral", "rnaviral_mode")]
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
    if "set_of_hmms" in cfg.__dict__:
            command.append(os.path.join(configs_dir, "hmm_mode.info"))


def prepare_config_spades(filename, cfg, log, additional_contigs_fname, K, stage, saves_dir, last_one, execution_home):
    subst_dict = dict()
    subst_dict["K"] = str(K)
    subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset)
    subst_dict["output_base"] = process_cfg.process_spaces(cfg.output_dir)
    subst_dict["tmp_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
    if additional_contigs_fname:
        subst_dict["additional_contigs"] = process_cfg.process_spaces(additional_contigs_fname)
        subst_dict["use_additional_contigs"] = bool_to_str(True)
    else:
        subst_dict["use_additional_contigs"] = bool_to_str(False)
    subst_dict["main_iteration"] = bool_to_str(last_one)
    subst_dict["entry_point"] = stage
    subst_dict["load_from"] = saves_dir
    if "checkpoints" in cfg.__dict__:
        subst_dict["checkpoints"] = cfg.checkpoints
    subst_dict["developer_mode"] = bool_to_str(cfg.developer_mode)
    subst_dict["time_tracer_enabled"] = bool_to_str(cfg.time_tracer)
    subst_dict["gap_closer_enable"] = bool_to_str(last_one or K >= options_storage.GAP_CLOSER_ENABLE_MIN_K)
    subst_dict["rr_enable"] = bool_to_str(last_one and cfg.rr_enable)
#    subst_dict["topology_simplif_enabled"] = bool_to_str(last_one)
    subst_dict["max_threads"] = cfg.max_threads
    subst_dict["max_memory"] = cfg.max_memory
    subst_dict["save_gp"] = bool_to_str(cfg.save_gp)
    if not last_one:
        subst_dict["correct_mismatches"] = bool_to_str(False)
    if "resolving_mode" in cfg.__dict__:
        subst_dict["resolving_mode"] = cfg.resolving_mode
    if "pacbio_mode" in cfg.__dict__:
        subst_dict["pacbio_test_on"] = bool_to_str(cfg.pacbio_mode)
        subst_dict["pacbio_reads"] = process_cfg.process_spaces(cfg.pacbio_reads)
    if cfg.cov_cutoff == "off":
        subst_dict["use_coverage_threshold"] = bool_to_str(False)
    else:
        subst_dict["use_coverage_threshold"] = bool_to_str(True)
        if cfg.cov_cutoff == "auto":
            subst_dict["coverage_threshold"] = 0.0
        else:
            subst_dict["coverage_threshold"] = cfg.cov_cutoff
    if cfg.lcer_cutoff is not None:
        subst_dict["lcer_enabled"] = bool_to_str(True)
        subst_dict["lcer_coverage_threshold"] = cfg.lcer_cutoff

    if "series_analysis" in cfg.__dict__:
        subst_dict["series_analysis"] = cfg.series_analysis
    process_cfg.substitute_params(filename, subst_dict, log)


def prepare_config_rnaspades(filename, log):
    if not options_storage.args.rna:
        return
    subst_dict = dict()
    subst_dict["ss_enabled"] = bool_to_str(options_storage.args.strand_specificity is not None)
    subst_dict["antisense"] = bool_to_str(options_storage.args.strand_specificity == "rf")
    process_cfg.substitute_params(filename, subst_dict, log)

def prepare_config_bgcspades(filename, cfg, log):
    if not "set_of_hmms" in cfg.__dict__:
        return
    subst_dict = dict()
    subst_dict["set_of_hmms"] = cfg.set_of_hmms
    if options_storage.args.bio:
        subst_dict["component_size_part"] = 1
        subst_dict["set_copynumber"] = bool_to_str(True)
        subst_dict["start_only_from_tips"] = bool_to_str(True)
    process_cfg.substitute_params(filename, subst_dict, log)

def prepare_config_construction(filename, log):
    if options_storage.args.read_cov_threshold is None:
        return
    subst_dict = dict()
    subst_dict["read_cov_threshold"] = options_storage.args.read_cov_threshold
    process_cfg.substitute_params(filename, subst_dict, log)


class IterationStage(stage.Stage):
    def __init__(self, K, prev_K, last_one, get_stage, latest, *args):
        super(IterationStage, self).__init__(*args)
        self.K = K
        self.short_name = "k%d" % self.K
        self.prev_K = prev_K
        self.last_one = last_one
        self.get_stage = get_stage
        self.latest = latest

    def generate_config(self, cfg):
        data_dir = os.path.join(cfg.output_dir, "K%d" % self.K)
        saves_dir = os.path.join(data_dir, "saves")
        dst_configs = os.path.join(data_dir, "configs")

        if self.get_stage(self.short_name) == options_storage.BASE_STAGE:
            if not os.path.isdir(data_dir):
                os.makedirs(data_dir)
            support.copy_tree(os.path.join(self.tmp_configs_dir, "debruijn"), dst_configs, preserve_times=False)

        if self.prev_K:
            additional_contigs_dname = os.path.join(cfg.output_dir, "K%d" % self.prev_K, "simplified_contigs")
        else:
            additional_contigs_dname = None

        if "read_buffer_size" in cfg.__dict__:
            # FIXME why here???
            process_cfg.substitute_params(os.path.join(dst_configs, "construction.info"),
                                          {"read_buffer_size": cfg.read_buffer_size}, self.log)
        if "scaffolding_mode" in cfg.__dict__:
            # FIXME why here???
            process_cfg.substitute_params(os.path.join(dst_configs, "pe_params.info"),
                                          {"scaffolding_mode": cfg.scaffolding_mode}, self.log)

        prepare_config_rnaspades(os.path.join(dst_configs, "rna_mode.info"), self.log)
        prepare_config_bgcspades(os.path.join(dst_configs, "hmm_mode.info"), cfg, self.log)
        prepare_config_construction(os.path.join(dst_configs, "construction.info"), self.log)
        cfg_fn = os.path.join(dst_configs, "config.info")
        prepare_config_spades(cfg_fn, cfg, self.log, additional_contigs_dname, self.K, self.get_stage(self.short_name),
                              saves_dir, self.last_one, self.bin_home)

    def get_command(self, cfg):
        data_dir = os.path.join(cfg.output_dir, "K%d" % self.K)
        dst_configs = os.path.join(data_dir, "configs")
        cfg_fn = os.path.join(dst_configs, "config.info")
        args = [cfg_fn]
        add_configs(args, dst_configs, cfg)

        command = [commands_parser.Command(
            STAGE="K%d" % self.K,
            path=os.path.join(self.bin_home, "spades-core"),
            args=args,
            config_dir=os.path.relpath(data_dir, options_storage.args.output_dir),
            short_name=self.short_name)]
        return command
