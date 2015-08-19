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
sys.path.append(os.path.join(sys.path[0], "../../../ext/src/python_libs/"))
from process_cfg import *
if sys.version.startswith('2.'):
    import pyyaml2 as pyyaml
elif sys.version.startswith('3.'):
    import pyyaml3 as pyyaml
import support

########################################################################

if len(sys.argv) != 2:
	print ("Converts datasets from .info format into new .yaml format (for each .info file output is modified .info file and new .yaml file)\n")
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
    yaml_file = os.path.join(working_dir, os.path.splitext(os.path.basename(info_file))[0] + ".yaml")  
    if os.path.isfile(yaml_file):
        print "Skipping", info_file, "because corresponding .yaml file already exists"
        continue
    print "\tProcessing", info_file,
    dataset_data = [{}, {}]
    content = load_dataset_from_info_file(info_file)
    if "reads" in content:
        print "\nSkipping", info_file, "because it contains link to a .yaml file! (" + content["reads"] + ")"
        continue
    for k, v in content.items():        
        if k.find("_reads") != -1 or k.find("jumping_") != -1:
            if k.find("jumping_") != -1:
                lib_type = 'mp'
            else:
                lib_type = 'pe'                
            if v.startswith('"'):
                v = v[1:-1].strip()
            reads = v.split()
            if k.startswith("paired_reads"):
                if len(reads) == 1:
                    add_to_dataset('--12', reads[0], dataset_data, lib_type)
                else:
                    add_to_dataset('-1', reads[0], dataset_data, lib_type)
                    add_to_dataset('-2', reads[1], dataset_data, lib_type)
            elif k.find("single") != -1:
                for read in reads:
                    add_to_dataset('-s', read, dataset_data, lib_type)
            elif k == "jumping_first":
                    add_to_dataset('-1', reads[0], dataset_data, lib_type)
            elif k == "jumping_second":
                    add_to_dataset('-2', reads[0], dataset_data, lib_type)
            else:
                print >> sys.stderr, "\nError: reads are not paired and not single!"
                continue

    dataset_data = support.correct_dataset(dataset_data)
    print '...writing to .yaml:', yaml_file, "and updating .info file with link to .yaml:", info_file
    pyyaml.dump(dataset_data, file(yaml_file, 'w'))
    info = open(info_file, 'a')
    info.write("\n")
    info.write("reads\t" + os.path.basename(yaml_file) + "\n")
    #print "yaml"
    #print pyyaml.dump(dataset_data)

