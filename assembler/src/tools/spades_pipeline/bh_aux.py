#!/usr/bin/env python

import os
import sys
import re

from process_cfg import bool_to_str

def verify(expr, message):
    if (not (expr)):
        print "Assertion failed. Message: " + message
        exit(0)

def check_files(prefix):
    msg = "Failure! Check output files."
    verify(os.path.isfile(prefix + ".cor.fastq"), msg + " prefix is " + prefix + " error 1")
    verify(os.path.isfile(prefix + ".unp.fastq"), msg + " prefix is " + prefix + " error 2")
    verify(os.path.isfile(prefix + ".bad.fastq"), msg + " prefix is " + prefix + " error 3")

def determine_it_count(tmp_dir):
    import re
    files = os.listdir(tmp_dir)
    answer = 0;
    for f in files:
        m = re.match(r"^(\d+)\..*", f)
        if m:
            val = int(m.group(1))
            if (val > answer):
                answer = val
    return answer            

def determine_read_files(folder, str_it_count, input_files):
    id = 1    
    answer = dict()
    answer["paired_reads"] = '"'
    answer["single_reads"] = '"'

    for input_file in input_files:
        prefix = os.path.basename(input_file) + "." + str_it_count        
        check_files(folder + prefix)

        answer["paired_reads"] += folder + prefix + ".cor.fastq" + '  '
        answer["single_reads"] += folder + prefix + ".unp.fastq" + '  '

    answer["paired_reads"] += '"'
    answer["single_reads"] += '"'
    
    return answer

# based on hammer function in ./src/tools/datasets.py
def hammer(given_props, output_dir, compress):
    def read_files():
        read_files = ["paired_reads", "single_reads"]        
        return read_files

    cmd = []
    dataset_entry = []
    for prop in given_props.iterkeys():        
        val = given_props[prop]
        if prop in read_files():
            new_val = '"'
            for oldfile in val[1:-1].strip().split("  "):
                newfile = output_dir + "/" + os.path.basename(oldfile)
                cmd += ["cp " + oldfile + " " + output_dir + "/"]
                if compress:
                    cmd += ["gzip -9 " + newfile]
                    newfile += ".gz"
                new_val += os.path.relpath(newfile, output_dir) + '  '
            new_val += '"'
            val = new_val
        dataset_entry += [(prop, val)]
    #dataset_entry = map(lambda (a, b): a + "\t" + b, dataset_entry)
    #dataset_entry = reduce(lambda x, y: x + "\n" + y, dataset_entry)
    for c in cmd:
        os.system(c)
    return dataset_entry

def dataset_pretty_print(dataset):    
    canonical_order       = ["paired_reads", "single_reads", "jumping_first", "jumping_second", "jumping_single_first", "jumping_single_second", "RL", "IS", "delta", "jump_is", "jump_rl", "single_cell", "is_var", "reference_genome", "genes", "operons"]
    max_property_name_len = len(max(canonical_order, key=len))
    tabulation            = "    " 
    
    pretty = ""
    dataset_dict = dict(dataset)
    for prop in canonical_order:
        if dataset_dict.has_key(prop):
            pretty += prop.ljust(max_property_name_len) + tabulation + dataset_dict[prop] + "\n"
    return pretty

def generate_dataset(cfg, tmp_dir, input_files):    
    str_it_count = str(determine_it_count(tmp_dir))

    if (cfg.max_iterations <= 10):
        str_it_count = "0" + str_it_count
    dataset_cfg = determine_read_files(tmp_dir + r"/", str_it_count, input_files)
    dataset_cfg["single_cell"] = bool_to_str(cfg.single_cell)
    if cfg.__dict__.has_key("reference"):    
        dataset_cfg["reference_genome"] = os.path.abspath(os.path.expandvars(cfg.reference))
    if cfg.__dict__.has_key("genes"):    
        dataset_cfg["genes"]            = os.path.abspath(os.path.expandvars(cfg.genes))
    if cfg.__dict__.has_key("operons"):        
        dataset_cfg["operons"]          = os.path.abspath(os.path.expandvars(cfg.operons))
    return dataset_pretty_print(hammer(dataset_cfg, cfg.output_dir, cfg.gzip_output))

def split_paired_file(cfg):
    input_filename = cfg.input_reads[0]
    ext = os.path.splitext(input_filename)[1]

    input_file = file    
    out_basename = ""
    
    if ext == '.gz':
        import gzip
        input_file   = gzip.open(input_filename, 'r')
        ungzipped    = os.path.splitext(input_filename)[0]
        out_basename = os.path.splitext(os.path.basename(ungzipped))[0]   
    else:
        input_file   = open(input_filename, 'r')
        out_basename = os.path.splitext(os.path.basename(input_filename))[0]

    out_left_filename  = os.path.join(cfg.working_dir, out_basename + "_left.fastq") 
    out_right_filename = os.path.join(cfg.working_dir, out_basename + "_right.fastq")    
    
    print("== Splitting " + input_filename + " into left and right reads")

    out_left_file = open(out_left_filename, 'w')
    out_right_file = open(out_right_filename, 'w')
    for id, line in enumerate(input_file):
        if id % 8 < 4:
            out_left_file.write(line)
        else: 
            out_right_file.write(line)
    
    out_left_file.close()
    out_right_file.close()
    input_file.close()

    return [out_left_filename, out_right_filename]

