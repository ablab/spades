#!/usr/bin/python

import sys
import os
import shutil
import re
import getopt
import datetime
import subprocess
import glob

###################################################################

sys.path.append(os.path.join(os.path.abspath(sys.path[0]), 'conversion'))
sys.path.append(os.path.join(os.path.abspath(sys.path[0]), 'stat'))

bowtie_path  = os.path.join(os.path.abspath(sys.path[0]), '../../../ext/tools/bowtie-0.12.7')
bowtie_build = os.path.join(bowtie_path, "bowtie-build")
bowtie       = os.path.join(bowtie_path, "bowtie")

tmp_folder = "tmp"
output_dir = "results_" + datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
thread_num = 16
make_latest_symlink = True

###################################################################

long_options = "output-dir= thread-num=".split()
short_options = "o:t:"

def usage():
    print 'Estimation reads quality'
    print 'Usage:', sys.argv[0], ' [options described below] datasets info-file(s) (for the same reference!)'
    print ""
    print "Options with parameters:"
    print "-o\t--output-dir\tDirectory to store all result files"
    print "-t\t--thread-num\tMax number of threads (default is " + str(thread_num) + ")"
    
def check_file(f):
    if not os.path.isfile(f):
        print "Error - file not found:", f
        sys.exit(2)
    return f

try:
    options, datasets = getopt.gnu_getopt(sys.argv[1:], short_options, long_options)
except getopt.GetoptError, err:
    print str(err)
    print ""
    usage()
    sys.exit(1)

for opt, arg in options:
    if opt in ('-o', "--output-dir"):
        output_dir = arg
        make_latest_symlink = False  
    elif opt in ('-t', "--thread-num"):
        thread_num = int(arg)
        if thread_num < 1:
            thread_num = 1      
    else:
        raise ValueError

for d in datasets:
    check_file(d)    

if not datasets:
    usage()
    sys.exit(1)   

###################################################################

def get_full_path(dataset, rel_path):
    return os.path.abspath(os.path.join(os.path.dirname(dataset), rel_path))

def ungzip_if_needed(filename, output_folder, force_copy = False):
    file_basename, file_extension = os.path.splitext(filename)
    if file_extension == ".gz":
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        ungzipped_filename = os.path.join(output_folder, os.path.basename(file_basename))
        ungzipped_file = open(ungzipped_filename, 'w')
        subprocess.call(['gunzip', filename, '-c'], stdout=ungzipped_file)
        ungzipped_file.close()
        filename = ungzipped_filename    
    elif force_copy:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        shutil.copy(filename, output_folder)
        filename = os.path.join(output_folder, os.path.basename(filename))
    return filename

###################################################################

reference = ""

datasets_dict = dict()

print("Analyzing datasets")
for dataset in datasets:

    basename = os.path.splitext(os.path.basename(dataset))[0]
    cur_key = basename
    i = 1
    while datasets_dict.has_key(cur_key):
        cur_key = basename + "_" + str(i)

    skip_this = False
    cur_reads = []
    cur_reference = ""
    for line in open(dataset, 'r'):
        if line.startswith("reference_genome"):
            cur_reference = get_full_path(dataset, line.split()[1].strip())
            if not reference:
                reference = cur_reference
            elif reference != cur_reference:
                skip_this = True
                break
        elif line.startswith("paired_reads") or line.startswith("single_reads"):
            line = line.replace('"', '')
            reads = line.split()[1:]
            for read in reads:
                cur_reads.append(get_full_path(dataset, read))            
    
    if not cur_reference:
        print("  " + dataset + " was skipped because it contains no reference")    
        continue       
    if skip_this:
        print("  " + dataset + " was skipped because of different reference (comparison only between all datasets of one organism)")    
        continue
    if len(cur_reads) == 0:
        print("  " + dataset + " was skipped because it contains no reads")    
        continue
        
    datasets_dict[cur_key] = cur_reads   
    print("  " + dataset + " ==> " + cur_key)

if not reference:
    print("Can't continue estimation - no reference in all datasets")
    sys.exit(1)

###################################################################

tmp_folder = os.path.join(output_dir, tmp_folder)
if not os.path.exists(tmp_folder):
    os.makedirs(tmp_folder)

print("Unpacking data (if needed) to temporary folder (" + tmp_folder + ") and N-trimming")

print("  reference...")
reference = ungzip_if_needed(reference, tmp_folder)

# TODO fastA analysis (we should convert all in fasta if there is at least one file in fasta)
import trim_ns
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    ungzipped_reads = []
    for read in datasets_dict[dataset]:    
        copied_read = ungzip_if_needed(read, os.path.join(tmp_folder, dataset), True)
        trim_ns.trim_file(copied_read, copied_read)
        ungzipped_reads.append(copied_read)
    datasets_dict[dataset] = ungzipped_reads

# creating index
index_folder = os.path.join(tmp_folder, "index")
if not os.path.exists(index_folder):
    os.makedirs(index_folder)
index_name   = os.path.splitext(os.path.basename(reference))[0]
index        = os.path.join(index_folder, index_name)
print("Creating index " + index)
index_log = open(os.path.join(output_dir, "index.log"),'w')
index_err = open(os.path.join(output_dir, "index.err"),'w')
subprocess.call([bowtie_build, reference, index], stdout=index_log, stderr=index_err)
index_log.close()
index_err.close()

# bowtie-ing
print("Aligning")
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    align_log = open(os.path.join(output_dir, dataset + ".log"),'w')
    align_err = open(os.path.join(output_dir, dataset + ".err"),'w') 
    reads_string = reduce(lambda x, y: x + ',' + y, datasets_dict[dataset])   
    subprocess.call([bowtie, '-c', '-q', '--suppress', '6,7,8', index, '-p', str(thread_num), reads_string], stdout=align_log, stderr=align_err)
    align_log.close()
    align_err.close()    

# TODO other logic
for k,v in datasets_dict.iteritems():
    print "k=", k, "v=", v
    





