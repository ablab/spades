#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil

import commands_parser
import options_storage
from stages import stage
import process_cfg
from process_cfg import bool_to_str

# TODO copypast from iteraton stage
READS_TYPES_USED_IN_CONSTRUCTION = ["paired-end", "single", "hq-mate-pairs"]
READS_TYPES_USED_IN_RNA_SEQ = ["paired-end", "single", "trusted-contigs", "untrusted-contigs"]


def add_configs(command, configs_dir):
    # Order matters here!
    mode_config_mapping = [("single_cell", "mda_mode"),
                           ("meta", "meta_mode"),
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

