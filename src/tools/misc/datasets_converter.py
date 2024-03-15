#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys

sys.path.append(os.path.join(sys.path[0], "../../spades_pipeline/"))
from process_cfg import *

########################################################################

# for pretty-printing
canonical_order       = ["paired_reads", "single_reads", "jumping_first", "jumping_second", "jumping_single_first", "jumping_single_second", "RL", "IS", "delta", "jump_is", "jump_rl", "single_cell", "is_var", "reference_genome"]
max_property_name_len = len(max(canonical_order, key=len))
tabulation            = "    " 

########################################################################

if len(sys.argv) != 4:
	print ("Splits old-style datasets.info file into separate files per each dataset (new format)\n")
	print ("Usage: " + sys.argv[0] + " old_style_datasets.info folder_for_separate_datasets folder_with_reads (for relative paths in new format datasets)")
	exit(0)

output_folder = sys.argv[2]
reads_folder  = sys.argv[3]

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

old_datasets = load_config_from_info_file(sys.argv[1])
for (key, value) in old_datasets.iteritems():
    if key != "common" and key != "OPPA": # oops, hard-code (for datasets_archive)
        cur_dataset = open(os.path.join(output_folder, key + ".info"), 'w')
        cur_ds_dict = dict()
        paired_reads = '"'
        single_reads = '"'
        for (prop, value) in value.__dict__.iteritems():
            if prop in ["RL", "IS", "delta", "jump_is", "jump_rl", "is_var"]:
                cur_ds_dict[prop] = str(value)                
            elif prop == "single_cell":
                cur_ds_dict[prop] = bool_to_str(value)
            else:                
                value = os.path.relpath(os.path.join(reads_folder, value), output_folder)
                if prop in ["first", "second"]:
                    paired_reads += value + ' '
                elif prop in ["single_first", "single_second"]:
                    single_reads += value + ' '
                else:
                    cur_ds_dict[prop] = value                

        paired_reads += '"'
        single_reads += '"'
        if paired_reads != '""':
            cur_ds_dict["paired_reads"] = paired_reads
        if single_reads != '""':
            cur_ds_dict["single_reads"] = single_reads
        
        # pretty-printing
        for prop in canonical_order:
            if cur_ds_dict.has_key(prop):
                cur_dataset.write(prop.ljust(max_property_name_len) + tabulation + cur_ds_dict[prop] + "\n")

        cur_dataset.close()
