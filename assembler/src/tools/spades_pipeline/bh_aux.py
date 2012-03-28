#!/usr/bin/env python

import os
import sys
import subprocess
import re

from process_cfg import bool_to_str

def verify(expr, message):
    if (not (expr)):
        print "Assertion failed. Message: " + message
        exit(0)

def check_files(files, str_it_count):
    msg = "Failure! Check output files."
    verify(len(files) == 6, msg + "1")
    verify(str_it_count + ".reads.0.left.corrected.fastq" in files, msg + "2")
    verify(str_it_count + ".reads.0.right.corrected.fastq" in files, msg + "3")
    verify(str_it_count + ".reads.0.left.unpaired.fastq" in files, msg + "4")
    verify(str_it_count + ".reads.0.right.unpaired.fastq" in files, msg + "5")
    verify(str_it_count + ".reads.0.left.bad.fastq" in files, msg + "6")
    verify(str_it_count + ".reads.0.right.bad.fastq" in files, msg + "7")

def determine_it_count(tmp_dir):
    import re
    files = subprocess.check_output('ls -1 ' + tmp_dir, shell=True).strip().split('\n')
    answer = 0;
    for f in files:
        m = re.match(r"^(\d+)\..*", f)
        if m:
            val = int(m.group(1))
            if (val > answer):
                answer = val
    return answer            

def determine_read_files(folder, str_it_count):
    files = subprocess.check_output('ls -1 ' + folder + r"/" + str_it_count + "*.reads.* | xargs -n1 basename"\
    , shell=True).strip().split('\n')
    check_files(files, str_it_count)
    
    answer = dict()    
    answer["first"] = folder + str_it_count + ".reads.0.left.corrected.fastq"    
    answer["second"] = folder + str_it_count + ".reads.0.right.corrected.fastq"
    answer["single_first"] = folder + str_it_count + ".reads.0.left.unpaired.fastq"    
    answer["single_second"] = folder + str_it_count + ".reads.0.right.unpaired.fastq"
    return answer

# based on hammer function in ./src/tools/datasets.py
def hammer(given_props, output_dir, compress):
    def read_files():
        read_files = ["first", "second"]
        read_files += ["single_" + x for x in read_files]
        return read_files

    cmd = []
    dataset_entry = []
    for prop in given_props.iterkeys():        
        val = given_props[prop]
        if prop in read_files():
            newfile = output_dir + "/" + os.path.basename(val)
            cmd += ["cp " + val + " " + output_dir + "/"]
            if compress:
                cmd += ["gzip -9 " + newfile]
                newfile += ".gz"
            val = os.path.relpath(newfile, output_dir) # re.sub(".*input/", "", newfile)
        dataset_entry += [(prop, val)]
    dataset_entry = map(lambda (a, b): a + "\t" + b, dataset_entry)
    dataset_entry = reduce(lambda x, y: x + "\n" + y, dataset_entry)
    for c in cmd:
        os.system(c)
    return dataset_entry

def generate_dataset(cfg, tmp_dir):    
    str_it_count = str(determine_it_count(tmp_dir))

    if (cfg.max_iterations <= 10):
        str_it_count = "0" + str_it_count
    dataset_cfg = determine_read_files(tmp_dir + r"/", str_it_count)
    dataset_cfg["single_cell"] = bool_to_str(cfg.single_cell)
    if cfg.__dict__.has_key("reference"):    
        dataset_cfg["reference_genome"] = os.path.abspath(os.path.expandvars(cfg.reference))
    if cfg.__dict__.has_key("genes"):    
        dataset_cfg["genes"]            = os.path.abspath(os.path.expandvars(cfg.genes))
    if cfg.__dict__.has_key("operons"):        
        dataset_cfg["operons"]          = os.path.abspath(os.path.expandvars(cfg.operons))
    return hammer(dataset_cfg, cfg.output_dir, cfg.gzip_output)
