#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import glob

sys.path.append(os.path.join(sys.path[0], "../../spades_pipeline/"))
from process_cfg import *
import support

########################################################################

if len(sys.argv) != 2:
	print ("Cleans datasets .info files from superfluous information\n")
	print ("Usage: " + sys.argv[0] + " <a DIRECTORY with datasets in .info format or a SINGLE .info FILE>")
	exit(0)

input_dir_or_file = sys.argv[1] 
if os.path.isdir(input_dir_or_file):
    working_dir = input_dir_or_file
    info_files = glob.glob(os.path.join(working_dir, "*.info"))
    if len(info_files) == 0:
        print (".info files not found in " + working_dir + " directory!")
        exit(0)
elif os.path.isfile(input_dir_or_file):
    info_files = [input_dir_or_file]
    working_dir = os.path.dirname(input_dir_or_file)
else:
    print (input_dir_or_file + " is not exist!")
    exit(0)

# aux function
def add_to_dataset(option, data, dataset_data, lib_type='pe'):
    data_type = support.get_data_type(option)    
    if lib_type == 'pe':
        record_id = 0
    else: # mate-pairs
        record_id = 1

    if not dataset_data[record_id]: # setting default values for a new record
        if lib_type == 'pe':
            dataset_data[record_id]['type'] = 'paired-end'
        else:
            dataset_data[record_id]['type'] = 'mate-pairs'            
    if data_type.endswith('reads'): # reads are stored as lists
        if data_type in dataset_data[record_id]:
            dataset_data[record_id][data_type].append(data)
        else:
            dataset_data[record_id][data_type] = [data]
    else: # other values are stored as plain strings
        dataset_data[record_id][data_type] = data


def load_dataset_from_info_file(info_file):
    content = dict()
    for line in open(info_file):
        if line.strip():
            key = line.split()[0]
            value = line[len(key) + 1:].strip()
            content[key] = value
    return content


for info_file in info_files:
    content = load_dataset_from_info_file(info_file)
    if "reads" in content:
        print '...updating .info file', info_file
        new_info = open(info_file, 'w')    
        new_info.write("reads\t" + content["reads"] + "\n")
        if "single_cell" in content:
            new_info.write("single_cell\t" + content["single_cell"] + "\n")
        if "reference_genome" in content:
            new_info.write("reference_genome\t" + content["reference_genome"] + "\n")
        new_info.write("\n")
        for info_field in ["RL", "IS", "delta", "jump_rl", "jump_is", "jump_delta"]:
            if info_field in content:
                new_info.write("; " + info_field + "\t" + content[info_field] + "\n")
        new_info.close()
    else:
        print "\nSkipping", info_file, "because it doesn't contains link to a .yaml file!\n"

